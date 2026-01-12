#!/usr/bin/env python3
from typing import Any, FrozenSet, Callable
import data_organizer
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

Key = tuple[str, FrozenSet[tuple[str, Any]]]

def make_key(name: str, params: dict[str, Any]) -> Key:
    return (name, frozenset(params.items()))

def get_key_name(key: Key) -> str:
    return key[0]

def get_key_params(key: Key) -> dict[str, Any]:
    return dict(key[1])

def exponential_ansatz(t, C, V):
    return C * np.exp(-V * t)

def cornell_potential_ansatz(r, A, sigma, B):
    return A + sigma * r - B / r

def register(name: str):
    def decorator(method):
        method._calc_name = name
        return method
    return decorator


class Calculator:
    _registry: dict[str, Callable] = {}

    @classmethod
    def _populate_registry(cls):
        """Scans the class for registered methods and populates the registry."""
        # Using dir(cls) is safer than __dict__ as it covers the MRO if you ever inherit
        for attr_name in dir(cls):
            attr = getattr(cls, attr_name)
            if callable(attr) and hasattr(attr, "_calc_name"):
                cls._registry[attr._calc_name] = attr

    def __init__(self, file_data: data_organizer.FileData, n_bootstrap: int = 200, seed: int = 42, step_size: int = 1):
        # Trigger registration if it hasn't happened yet
        if not self._registry:
            self._populate_registry()

        self.file_data = file_data
        self.variables: dict[Key, data_organizer.VariableData] = {}
        
        # Bootstrap configuration
        self.n_bootstrap = n_bootstrap
        self.seed = seed
        self.step_size = step_size
        
        # Cache for indices
        self._bootstrap_indices: np.ndarray | None = None

    def get_variable(self, name: str, **params) -> data_organizer.VariableData:
        key = make_key(name, params)
        if key in self.variables:
            return self.variables[key]
        if name in self._registry:
            try:
                # Call the registered method (passing self and params)
                var_data = self._registry[name](self, **params)
            except TypeError as e:
                raise TypeError(f"Error calculating '{name}'. Check if params {params} match method signature.") from e
            
            self.variables[key] = var_data
            return var_data

        raise KeyError(f"No calculator registered for variable '{name}'")
    
    def get_observable(self, obs_name: str) -> data_organizer.ObservableData:
        for obs in self.file_data.observables:
            if obs.name == obs_name:
                return obs
        raise KeyError(f"Observable '{obs_name}' not found in file data.")
    
    def _get_bootstrap_indices(self, n_samples: int) -> np.ndarray:
        """
        Returns (n_bootstrap, n_effective) indices.
        Reduces the sample pool by step_size to mitigate autocorrelation.
        """
        n_effective = n_samples // self.step_size
        
        if self._bootstrap_indices is None or self._bootstrap_indices.shape[1] != n_effective:
            rng = np.random.default_rng(self.seed)
            reduced_indices = rng.integers(0, n_effective, size=(self.n_bootstrap, n_effective))
            self._bootstrap_indices = reduced_indices * self.step_size
        
        return self._bootstrap_indices
    
    ### Variable implementations ###

    @register("W_R_T")
    def _calc_W_R_T(self, R: int, T: int) -> data_organizer.VariableData:
        try:
            obs_val = self.file_data.get("W_temp").values 
            obs_L   = self.file_data.get("L").values
            obs_T   = self.file_data.get("T").values
        except ValueError as e:
            raise KeyError(f"Missing required columns (L, T, or W_temp) in file: {e}")

        arr_val = np.array(obs_val)
        arr_L   = np.array(obs_L)
        arr_T   = np.array(obs_T)

        # 3. Create a Mask (Where L==R AND T==T)
        mask = (arr_L == R) & (arr_T == T)
        selected_values = arr_val[mask]
        n_samples = len(selected_values)

        if n_samples == 0:
            raise ValueError(f"No data found for R={R}, T={T}")

        var_data = data_organizer.VariableData("W_R_T")
        
        mean_val = np.mean(selected_values)

        # 3. Bootstrap Calculation
        indices = self._get_bootstrap_indices(n_samples) # Shape: (n_boot, n_samples)
        bootstrapped_rows = selected_values[indices] 
        bootstrap_means = np.mean(bootstrapped_rows, axis=1)

        var_data.set_value(mean_val, bootstrap_samples=bootstrap_means, R=R, T=T)
        return var_data
    
    @register("V_R")
    def _calc_V_R(self, R: int, t_min: int = 1) -> data_organizer.VariableData:
        """
        Calculates V(R) and its error by fitting W(R,T).
        Includes robust handling for negative signal and array shape mismatches.
        """
        all_L = np.array(self.file_data.get("L").values)
        all_T = np.array(self.file_data.get("T").values)
        unique_ts = np.unique(all_T[all_L == R])

        # Filter T >= t_min
        ts_to_fit = sorted([t for t in unique_ts if t >= t_min])
        
        # W_R_T Variables
        w_vars = [self.get_variable("W_R_T", R=R, T=t) for t in ts_to_fit]
        
        ws_main = np.array([v.get() for v in w_vars])
        ts_fit  = np.array(ts_to_fit)
        ws_err  = np.array([v.err() for v in w_vars])

        # --- Positivity Cut ---
        # Find the first index where data becomes non-positive
        bad_indices = np.where(ws_main <= 0)[0]
        if len(bad_indices) > 0:
            first_bad = bad_indices[0]
            # Truncate everything after and including the first bad point
            ws_main = ws_main[:first_bad]
            ts_fit  = ts_fit[:first_bad]
            ws_err  = ws_err[:first_bad]
            
        # Check if we still have enough points after truncation
        if len(ts_fit) < 3:
            # Signal lost (values became negative) too early.
            var_data = data_organizer.VariableData("V_R")
            var_data.set_value(np.nan, bootstrap_samples=[np.nan]*self.n_bootstrap, R=R)
            return var_data

        # Store the valid length to slice bootstrap samples later
        valid_len = len(ts_fit)

        # Handle case where errors might be zero
        if np.any(ws_err == 0):
             mean_err = np.mean(ws_err[ws_err > 0]) if np.any(ws_err > 0) else 1.0
             ws_err[ws_err == 0] = mean_err

        def fit_potential(ts, ws, sigma=None):
            # Check for empty or non-positive data before fitting to avoid log errors
            if len(ws) == 0 or np.any(ws <= 0):
                return np.nan

            # initial guess
            try:
                denom = ts[1] - ts[0]
                if denom == 0: denom = 1.0
                p0_V = -np.log(ws[1]/ws[0]) / denom
                p0_C = ws[0] * np.exp(p0_V * ts[0])
            except (IndexError, ValueError, RuntimeWarning):
                p0_V, p0_C = 0.5, 1.0
                
            try:
                popt, _ = curve_fit(
                    exponential_ansatz, ts, ws, 
                    p0=[p0_C, p0_V], 
                    bounds=([-np.inf, -10.0], [np.inf, np.inf]),
                    sigma=sigma,
                    absolute_sigma=(sigma is not None)
                )
                return popt[1] # Return V
            except (RuntimeError, ValueError, TypeError):
                return np.nan

        # Weighted Fit
        V_main = fit_potential(ts_fit, ws_main, sigma=ws_err)
        
        if np.isnan(V_main):
            # If the main fit fails, return NaN immediately
            var_data = data_organizer.VariableData("V_R")
            var_data.set_value(np.nan, bootstrap_samples=[np.nan]*self.n_bootstrap, R=R)
            return var_data

        # Bootstrap Fits
        n_boot = self.n_bootstrap
        bootstrap_Vs = []

        all_boots = np.array([v.bootstrap() for v in w_vars])
        all_boots_T = all_boots.T 
        
        for i in range(n_boot):
            # FIX: Slice the bootstrap sample to match the truncated time array
            ws_sample = all_boots_T[i][:valid_len] 
            
            # Skip this universe if it fluctuates into negative values within the window
            if np.any(ws_sample <= 0):
                bootstrap_Vs.append(np.nan)
                continue
            
            # Fit using the SAME weights (ws_err) from the main fit
            v_boot = fit_potential(ts_fit, ws_sample, sigma=ws_err)
            
            # Append result (NaN or float)
            bootstrap_Vs.append(v_boot)

        var_data = data_organizer.VariableData("V_R")
        var_data.set_value(V_main, bootstrap_samples=bootstrap_Vs, R=R, t_min=t_min)
        return var_data
    
    @register("tau_int")
    def _calc_tau_int(self, obs_name: str = "plaquette", S: float = 1.5) -> data_organizer.VariableData:
        try:
            obs = self.get_observable(obs_name)
        except KeyError:
            # Fallback for commonly used names if exact match fails
            found = False
            for alias in ["plaquette", "retrace", "W_temp"]:
                 if alias in self.file_data.observables: # Simplified check
                     for o in self.file_data.observables:
                         if o.name == alias:
                             obs = o
                             found = True
                             break
            if not found:
                 raise KeyError(f"Observable '{obs_name}' not found for tau_int calculation.")

        series = np.array(obs.values)
        n = len(series)
        
        # Default fallback if series is too short
        tau = 0.5 

        if n >= 100:
            # Subtract mean
            series_centered = series - np.mean(series)
            var = np.var(series_centered)
            
            if var > 0:
                # FFT for autocorrelation
                ft = np.fft.rfft(series_centered)
                spec = np.abs(ft)**2
                acf = np.fft.irfft(spec)
                acf = acf[:n//2]
                
                # Normalize (N-normalization for FFT estimate)
                rho = (acf / np.arange(n, n - len(acf), -1)) / var

                # Integrate rho to find tau with automatic windowing
                tau_sum = 0.5
                for t in range(1, len(rho)):
                    tau_sum += rho[t]
                    if t >= S * tau_sum:
                        tau = tau_sum
                        break

        var_data = data_organizer.VariableData("tau_int")
        var_data.set_value(tau, obs_name=obs_name)
        return var_data

    @register("r0")
    def _calc_r0(self, t_min: int = 2, target_force: float = 1.65) -> data_organizer.VariableData:
        """
        Calculates the Sommer parameter r0/a by fitting the Cornell potential to V(R).
        Replaces analyze_wilson.fit_sommer_parameter.
        """
        # 1. Identify available R values
        # We assume R corresponds to 'L' in W_temp data
        try:
            all_L = np.array(self.file_data.get("L").values)
            unique_Rs = sorted(np.unique(all_L))
        except ValueError:
             raise KeyError("Could not determine available R (L) values from file data.")

        # 2. Gather V(R) for all R
        rs = []
        vs = []
        errs = []
        
        # We need a consistent set of bootstrap samples to propagate errors correctly
        # So we collect the bootstrap arrays from each V_R
        bootstrap_matrix = []

        for r in unique_Rs:
            try:
                # Reuse the existing V_R calculator
                v_var = self.get_variable("V_R", R=int(r), t_min=t_min)
                
                val = v_var.get()
                if np.isnan(val) or np.isinf(val): continue

                rs.append(r)
                vs.append(val)
                errs.append(v_var.err() if v_var.err() is not None else 1.0)
                
                # Collect bootstraps. Must ensure they align (same seed/length).
                # The Calculator class ensures consistent seeding in __init__.
                if v_var.bootstrap() is not None:
                    bootstrap_matrix.append(v_var.bootstrap())
                
            except (ValueError, RuntimeError):
                continue
        
        if len(rs) < 3:
            raise ValueError(f"Not enough valid V(R) points to fit r0 (found {len(rs)})")

        rs = np.array(rs)
        vs = np.array(vs)
        errs = np.array(errs)
        
        # Handle zero errors for fitting
        sigma = errs.copy()
        if np.any(sigma <= 0):
            mean_err = np.mean(sigma[sigma > 0]) if np.any(sigma > 0) else 1.0
            sigma[sigma <= 0] = mean_err

        def perform_fit(r_vals, v_vals, sigma_vals=None):
            try:
                popt, _ = curve_fit(
                    cornell_potential_ansatz, r_vals, v_vals, 
                    p0=[0.0, 0.1, 0.26], sigma=sigma_vals, 
                    absolute_sigma=(sigma_vals is not None), maxfev=2000
                )
                A, sig, B = popt
                
                # r0 definition: r^2 * (sigma + B/r^2) = 1.65  => r^2 = (1.65 - B)/sigma
                numerator = target_force - B
                if numerator < 0 or sig <= 0:
                    return np.nan, (A, sig, B)
                return np.sqrt(numerator / sig), (A, sig, B)
            except (RuntimeError, ValueError, TypeError):
                return np.nan, (0,0,0)

        # Main Fit
        r0_val, fit_params = perform_fit(rs, vs, sigma)
        corn_params = {'A': fit_params[0], 'sigma': fit_params[1], 'B': fit_params[2]}

        # Bootstrap Fits
        r0_bootstraps = []
        if len(bootstrap_matrix) == len(rs) and len(bootstrap_matrix) > 0:
            boot_data = np.array(bootstrap_matrix).T
            for i in range(len(boot_data)):
                val, _ = perform_fit(rs, boot_data[i], sigma)
                if not np.isnan(val):
                    r0_bootstraps.append(val)

        var_data = data_organizer.VariableData("r0")
        var_data.set_value(r0_val, bootstrap_samples=r0_bootstraps, t_min=t_min, cornell_params=corn_params)
        return var_data

    # --- Implementation of Creutz Ratio Analysis ---
    # The following methods implement the calculation of the lattice scale 'a'
    # using Creutz Ratios, as discussed in CreutzRatio.md (arXiv:2410.02794v1).
    #
    # Regarding the question "Is it possible to do it without smearing?":
    # Yes, it is possible. The calculations can be performed on unsmeared gauge
    # configurations (which corresponds to zero smearing/flow time). However, as the
    # paper points out, smearing significantly reduces the statistical noise (variance)
    # in the Wilson loops and Creutz ratios. Without smearing, you would need
    # much higher statistics (more measurements) to achieve a similar level of
    # precision for the calculated force and the resulting lattice scale r0.

    @register("chi")
    def _calc_chi(self, R: float, T: float) -> data_organizer.VariableData:
        """
        Calculates the Creutz ratio chi(R, T) using Wilson loops.
        R and T are the half-integer coordinates of the center of the Creutz plaquette.
        e.g., R=1.5, T=1.5 for the ratio of W(1,1), W(2,1), W(1,2), W(2,2).
        """
        # Creutz ratio is defined at half-integer coordinates
        # chi(t + a/2, r + a/2), so r_int = r - 0.5, t_int = t - 0.5
        r_int = int(R - 0.5)
        t_int = int(T - 0.5)
        
        try:
            # Get the four Wilson loops needed for the ratio
            w_tr   = self.get_variable("W_R_T", R=r_int,     T=t_int)
            w_t1_r = self.get_variable("W_R_T", R=r_int,     T=t_int + 1)
            w_t_r1 = self.get_variable("W_R_T", R=r_int + 1, T=t_int)
            w_t1_r1= self.get_variable("W_R_T", R=r_int + 1, T=t_int + 1)
        except (ValueError, KeyError) as e:
            raise ValueError(f"Could not get required Wilson loops for chi({R}, {T}): {e}")

        def calculate_log_ratio(w1, w2, w3, w4):
            # w1 = W(t+a, r), w2 = W(t, r+a), w3 = W(t,r), w4 = W(t+a, r+a)
            # Formula is ln( (w1 * w2) / (w3 * w4) ) assuming a=1
            num = w1 * w2
            den = w3 * w4
            if num > 0 and den > 0:
                return np.log(num / den)
            return np.nan

        # Main value
        chi_val = calculate_log_ratio(w_t1_r.get(), w_t_r1.get(), w_tr.get(), w_t1_r1.get())
        
        # Bootstrap values
        chi_boots = []
        # Get bootstrap samples for all 4 loops
        b_t1_r = w_t1_r.bootstrap()
        b_t_r1 = w_t_r1.bootstrap()
        b_tr   = w_tr.bootstrap()
        b_t1_r1= w_t1_r1.bootstrap()

        # Check if bootstrap samples are available
        if b_t1_r is None or b_t_r1 is None or b_tr is None or b_t1_r1 is None:
             raise ValueError(f"Bootstrap samples not available for one of the Wilson loops for chi({R},{T})")

        for i in range(self.n_bootstrap):
            boot_val = calculate_log_ratio(b_t1_r[i], b_t_r1[i], b_tr[i], b_t1_r1[i])
            chi_boots.append(boot_val)
        
        var_data = data_organizer.VariableData("chi")
        var_data.set_value(chi_val, bootstrap_samples=chi_boots, R=R, T=T)
        return var_data

    @register("F_chi")
    def _calc_F_chi(self, R: float, t_large: int) -> data_organizer.VariableData:
        """
        Extracts the static force F(R) from Creutz ratios chi(R, T)
        by taking the value at a large time extent, T = t_large.
        R is a half-integer (e.g., 1.5, 2.5).
        t_large is an integer, so T for chi will be t_large + 0.5.
        """
        T = float(t_large) + 0.5
        try:
            chi_var = self.get_variable("chi", R=R, T=T)
            var_data = data_organizer.VariableData("F_chi")
            # The value and bootstrap samples are directly from chi(R, t_large)
            var_data.set_value(chi_var.get(), bootstrap_samples=chi_var.bootstrap(), R=R, t_large=t_large)
            return var_data
        except (ValueError, KeyError) as e:
            raise ValueError(f"Could not calculate F_chi for R={R} at t_large={t_large}: {e}")

    @register("r0_chi")
    def _calc_r0_chi(self, t_large: int = 4, target_force: float = 1.65, max_rel_err: float = 0.5, use_weighted_fit: bool = True, fit_window: int = 2) -> data_organizer.VariableData:
        """
        Calculates Sommer parameter r0/a from Creutz ratios.
        It interpolates r^2 * F(r) to find where it equals target_force.
        Points with a high relative error on F(r) can be filtered out.
        A weighted linear fit over a window of points can be used for interpolation.
        """
        # 1. Identify available R values from data.
        try:
            all_L = np.array(self.file_data.get("L").values)
            unique_Ls = sorted(np.unique(all_L))
            # Rs for chi are like 1.5, 2.5, ... from L=1,2,3...
            unique_Rs = [float(L) + 0.5 for L in unique_Ls if L+1 in unique_Ls]
        except ValueError:
             raise KeyError("Could not determine available R (L) values from file data for r0_chi.")

        # 2. Gather F_chi(R) for all available R at a fixed large T.
        rs, fs, errs, bootstrap_matrix = [], [], [], []
        
        for r_val in unique_Rs:
            try:
                f_var = self.get_variable("F_chi", R=r_val, t_large=t_large)
                val, err = f_var.get(), f_var.err()

                if np.isnan(val) or np.isinf(val) or np.isnan(err) or np.isinf(err):
                    continue
                
                # Filter points with large relative error
                if max_rel_err is not None and val != 0 and abs(err / val) > max_rel_err:
                    continue

                rs.append(r_val)
                fs.append(val)
                errs.append(err)
                if f_var.bootstrap() is not None:
                    bootstrap_matrix.append(f_var.bootstrap())
            except (ValueError, KeyError):
                continue

        if len(rs) < 2: # Need at least 2 points for interpolation/fit
            raise ValueError(f"Not enough valid F_chi(R) points to process for r0 (found {len(rs)}).")
        
        rs, fs, errs = np.array(rs), np.array(fs), np.array(errs)

        def get_r0_from_force(r_vals, f_vals, f_errs=None):
            """
            Finds r0 by linearly interpolating/fitting r^2*F(r).
            Finds the first crossover and can use a weighted fit in a window around it.
            """
            if len(r_vals) < 2: return np.nan
            
            r2_vals = r_vals**2
            r2f_vals = r2_vals * f_vals

            # Find the first crossover point, assuming r_vals are sorted.
            for i in range(len(r_vals) - 1):
                y1, y2 = r2f_vals[i], r2f_vals[i+1]

                if (y1 <= target_force <= y2) or (y2 <= target_force <= y1):
                    # Crossover is between index i and i+1
                    if not use_weighted_fit or f_errs is None or len(f_errs) != len(r_vals):
                        # Simple linear interpolation
                        if y2 == y1: continue
                        r1_sq, r2_sq = r2_vals[i], r2_vals[i+1]
                        r0_sq = r1_sq + (target_force - y1) * (r2_sq - r1_sq) / (y2 - y1)
                        if r0_sq > 0: return np.sqrt(r0_sq)
                    else:
                        # Weighted linear fit over a window of points.
                        # Window includes `fit_window-1` points before and `fit_window` points after `i`.
                        start = max(0, i - fit_window + 1)
                        end = min(len(r_vals), i + fit_window + 1)
                        
                        win_r2 = r2_vals[start:end]
                        win_r2f = r2f_vals[start:end]
                        win_f_errs = f_errs[start:end]

                        if len(win_r2) < 2: continue

                        # Weights are inverse of variance of r^2*F(r) -> 1 / (r^4 * err(F)^2)
                        win_r4f_err2 = (win_r2**2) * (win_f_errs**2)
                        valid_w = win_r4f_err2 > 0
                        if np.sum(valid_w) < 2: continue # Not enough points with valid weights

                        weights = np.zeros_like(win_r4f_err2)
                        weights[valid_w] = 1.0 / win_r4f_err2[valid_w]

                        try:
                            # y = p[0]*x + p[1]
                            p = np.polyfit(win_r2[valid_w], win_r2f[valid_w], 1, w=weights[valid_w])
                            slope, intercept = p[0], p[1]
                        except (np.linalg.LinAlgError, ValueError):
                            continue # Fit failed

                        if slope == 0: continue
                        # Solve for r0^2: slope * r0^2 + intercept = target_force
                        r0_sq = (target_force - intercept) / slope
                        if r0_sq > 0: return np.sqrt(r0_sq)
            
            return np.nan # No crossover found

        # Main value
        r0_val = get_r0_from_force(rs, fs, errs)

        # Bootstrap
        r0_bootstraps = []
        if len(bootstrap_matrix) > 0 and len(bootstrap_matrix) == len(rs):
            boot_data = np.array(bootstrap_matrix).T
            for sample_fs in boot_data:
                # For bootstrap samples, we don't have individual errors,
                # so we use the errors from the main data set as weights.
                val = get_r0_from_force(rs, sample_fs, errs)
                r0_bootstraps.append(val)
        
        var_data = data_organizer.VariableData("r0_chi")
        var_data.set_value(r0_val, bootstrap_samples=r0_bootstraps, t_large=t_large)
        return var_data