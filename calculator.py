#!/usr/bin/env python3
from typing import Any, FrozenSet, Callable
import data_organizer
import numpy as np
from scipy.optimize import curve_fit

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
        Calculates V(R) and its error by fitting W(R,T) for every bootstrap sample.
        Performs a weighted fit using the errors of W(R,T) as weights.
        """
        all_L = np.array(self.file_data.get("L").values)
        all_T = np.array(self.file_data.get("T").values)
        unique_ts = np.unique(all_T[all_L == R])
        
        # Filter T >= t_min
        ts_to_fit = sorted([t for t in unique_ts if t >= t_min])
        
        if len(ts_to_fit) < 3:
            raise ValueError(f"Not enough time points to fit V(R) for R={R} (found {len(ts_to_fit)})")

        # W_R_T Variables
        w_vars = [self.get_variable("W_R_T", R=R, T=t) for t in ts_to_fit]
        
        ws_main = np.array([v.get() for v in w_vars])
        ts_fit  = np.array(ts_to_fit)
        ws_err = np.array([v.err() for v in w_vars])

        # Handle case where errors might be zero (to prevent fit issues)
        if np.any(ws_err == 0):
             # Replace 0 with mean of non-zero errors, or 1.0 if all are 0
             mean_err = np.mean(ws_err[ws_err > 0]) if np.any(ws_err > 0) else 1.0
             ws_err[ws_err == 0] = mean_err

        def fit_potential(ts, ws, sigma=None):
            # initial guess
            if ws[0] > 0 and ws[1] > 0:
                denom = ts[1] - ts[0]
                if denom == 0: denom = 1.0
                p0_V = -np.log(ws[1]/ws[0]) / denom
                p0_C = ws[0] * np.exp(p0_V * ts[0])
            else:
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
                print(f"Fit failed for V_R (R={R}) with ts={ts} and ws={ws}")
                return np.nan

        # Weighted Fit
        V_main = fit_potential(ts_fit, ws_main, sigma=ws_err)
        
        if np.isnan(V_main):
            raise ValueError(f"Main fit failed for V_R (R={R})")

        # Bootstrap Fits
        n_boot = self.n_bootstrap
        bootstrap_Vs = []

        all_boots = np.array([v.bootstrap() for v in w_vars])
        all_boots_T = all_boots.T 
        
        for i in range(n_boot):
            ws_sample = all_boots_T[i] # The W(T) curve for the i-th bootstrap universe
            
            # Fit this universe using the SAME weights (ws_err) from the main fit
            v_boot = fit_potential(ts_fit, ws_sample, sigma=ws_err)
            
            if not np.isnan(v_boot):
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