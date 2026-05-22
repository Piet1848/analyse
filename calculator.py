#!/usr/bin/env python3
import collections
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


def _positive_log_ratio(numerator, denominator):
    numerator_arr, denominator_arr = np.broadcast_arrays(
        np.asarray(numerator, dtype=float),
        np.asarray(denominator, dtype=float),
    )
    valid = (
        np.isfinite(numerator_arr)
        & np.isfinite(denominator_arr)
        & (numerator_arr > 0)
        & (denominator_arr > 0)
    )

    if numerator_arr.ndim == 0:
        if bool(valid):
            return float(np.log(float(numerator_arr) / float(denominator_arr)))
        return np.nan

    result = np.full(numerator_arr.shape, np.nan, dtype=float)
    result[valid] = np.log(numerator_arr[valid] / denominator_arr[valid])
    return result


def creutz_chi_from_wilson(w_rt, w_rt1, w_r1t, w_r1t1):
    """
    Standard adjacent-loop Creutz ratio.

    Callers pass Wilson-loop means or bootstrap means, matching the
    example_analysis convention of first grouping W_temp by (flow_time, L, T).
    """
    return _positive_log_ratio(
        np.asarray(w_rt1, dtype=float) * np.asarray(w_r1t, dtype=float),
        np.asarray(w_rt, dtype=float) * np.asarray(w_r1t1, dtype=float),
    )


def fit_r0_from_potential_data(
    rs,
    vs,
    errs=None,
    bootstrap_matrix=None,
    target_force: float = 1.65,
    n_bootstrap: int | None = None,
):
    """
    Fit the Cornell potential to preselected V(R) points and derive r0/a.

    Parameters
    ----------
    rs, vs:
        One-dimensional arrays with matching lengths.
    errs:
        Optional per-point uncertainties for weighted fitting.
    bootstrap_matrix:
        Optional array with shape (n_points, n_bootstrap). Each column is one
        bootstrap replica of the selected potential points.
    target_force:
        Sommer target used in the definition r0^2 * F(r0) = target_force.
    n_bootstrap:
        Number of bootstrap replicas to return when bootstrap_matrix is not
        provided. Defaults to the second dimension of bootstrap_matrix or 0.
    """
    rs_arr = np.asarray(rs, dtype=float)
    vs_arr = np.asarray(vs, dtype=float)

    if rs_arr.ndim != 1 or vs_arr.ndim != 1 or len(rs_arr) != len(vs_arr):
        raise ValueError("rs and vs must be one-dimensional arrays with matching lengths")

    if errs is None:
        sigma_arr = np.ones(len(rs_arr), dtype=float)
    else:
        sigma_arr = np.asarray(errs, dtype=float)
        if sigma_arr.ndim != 1 or len(sigma_arr) != len(rs_arr):
            raise ValueError("errs must be a one-dimensional array with the same length as rs")

    boot_arr = None
    if bootstrap_matrix is not None:
        boot_arr = np.asarray(bootstrap_matrix, dtype=float)
        if boot_arr.ndim != 2 or boot_arr.shape[0] != len(rs_arr):
            raise ValueError(
                "bootstrap_matrix must have shape (n_points, n_bootstrap) matching rs"
            )
        if n_bootstrap is None:
            n_bootstrap = int(boot_arr.shape[1])

    if n_bootstrap is None:
        n_bootstrap = 0

    valid_mask = np.isfinite(rs_arr) & np.isfinite(vs_arr) & np.isfinite(sigma_arr)
    fit_rs = rs_arr[valid_mask]
    fit_vs = vs_arr[valid_mask]
    fit_sigma = sigma_arr[valid_mask]

    if boot_arr is not None:
        boot_arr = boot_arr[valid_mask, :]

    if len(fit_rs) < 3:
        raise ValueError(f"Not enough valid V(R) points to fit r0 (found {len(fit_rs)}).")

    if np.any(fit_sigma <= 0):
        mean_err = np.mean(fit_sigma[fit_sigma > 0]) if np.any(fit_sigma > 0) else 1.0
        fit_sigma = fit_sigma.copy()
        fit_sigma[fit_sigma <= 0] = mean_err

    def perform_fit(r_vals, v_vals, sigma_vals=None):
        try:
            sigma_guess = (
                (v_vals[-1] - v_vals[0]) / (r_vals[-1] - r_vals[0])
                if len(r_vals) > 1 and r_vals[-1] != r_vals[0]
                else 0.1
            )
            B_guess = 0.26
            A_guess = v_vals[0] - sigma_guess * r_vals[0] + (B_guess / r_vals[0])
            popt, _ = curve_fit(
                cornell_potential_ansatz,
                r_vals,
                v_vals,
                p0=[A_guess, max(1e-4, sigma_guess), B_guess],
                sigma=sigma_vals,
                absolute_sigma=(sigma_vals is not None),
                maxfev=5000,
            )
            A, sig, B = popt
            numerator = target_force - B
            if numerator < 0 or sig <= 0:
                return np.nan, (A, sig, B)
            return np.sqrt(numerator / sig), (A, sig, B)
        except (RuntimeError, ValueError, TypeError):
            return np.nan, (0.0, 0.0, 0.0)

    r0_val, fit_params = perform_fit(fit_rs, fit_vs, fit_sigma)
    cornell_params = {"A": fit_params[0], "sigma": fit_params[1], "B": fit_params[2]}
    fit_model = cornell_potential_ansatz(fit_rs, *fit_params)
    residuals = fit_vs - fit_model
    chi2 = float(np.sum((residuals / fit_sigma) ** 2)) if len(fit_sigma) > 0 else None
    dof = max(int(len(fit_rs)) - len(fit_params), 0)
    chi2_dof = float(chi2 / dof) if chi2 is not None and dof > 0 else None

    fill_value = float(r0_val) if np.isfinite(r0_val) else np.nan
    r0_bootstraps = np.full(int(n_bootstrap), fill_value, dtype=float)
    if boot_arr is not None and boot_arr.size > 0:
        for i in range(boot_arr.shape[1]):
            repaired_vs = np.array(boot_arr[:, i], copy=True)
            invalid_points = ~np.isfinite(repaired_vs)
            if np.any(invalid_points):
                repaired_vs[invalid_points] = fit_vs[invalid_points]

            val, _ = perform_fit(fit_rs, repaired_vs, fit_sigma)
            r0_bootstraps[i] = fill_value if np.isnan(val) else val

    finite_boots = r0_bootstraps[np.isfinite(r0_bootstraps)]
    r0_err = float(np.std(finite_boots)) if finite_boots.size > 0 else None

    return {
        "r0": float(r0_val) if np.isfinite(r0_val) else np.nan,
        "r0_err": r0_err,
        "bootstrap_samples": r0_bootstraps,
        "cornell_params": cornell_params,
        "fit_rs": fit_rs,
        "fit_vs": fit_vs,
        "fit_errs": fit_sigma,
        "chi2": chi2,
        "dof": dof,
        "chi2_dof": chi2_dof,
    }

def calculate_volume_from_r0(
    lattice_extent: int | float,
    r0_on_a: float,
    r0_phys: float = 0.5,
) -> float:
    if r0_on_a is None or np.isnan(r0_on_a) or r0_on_a <= 0:
        return np.nan
    a = r0_phys / r0_on_a
    return lattice_extent * (a ** 4)

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
        self.step_size = int(step_size)
        if self.step_size < 1:
            raise ValueError("step_size must be at least 1")

        # Wilson-loop lookup is built lazily. For compact combined data we can reuse
        # the pre-grouped arrays directly and avoid rebuilding the cache.
        self._w_rt_cache: dict[data_organizer.FlowKey, np.ndarray] | None = None
        self._pair_order: list[data_organizer.FlowKey] | None = None

        compact_wilson = getattr(self.file_data, "wilson_by_flow_pair", None)
        if compact_wilson is not None:
            self._w_rt_cache = dict(compact_wilson)
            pair_order = getattr(self.file_data, "flow_pair_order", None)
            if pair_order is not None:
                self._pair_order = [
                    (data_organizer._normalize_flow_time(flow_time), int(r), int(t))
                    for flow_time, r, t in pair_order
                ]

    def _canonicalize_params(self, params: dict[str, Any]) -> dict[str, Any]:
        canonical: dict[str, Any] = {}
        for name, value in params.items():
            if name == "flow_time":
                normalized_flow = data_organizer._normalize_flow_time(value)
                if normalized_flow is None:
                    continue
                canonical[name] = normalized_flow
            elif isinstance(value, np.integer):
                canonical[name] = int(value)
            elif isinstance(value, np.floating):
                canonical[name] = float(value)
            else:
                canonical[name] = value
        return canonical

    def get_variable(self, name: str, **params) -> data_organizer.VariableData:
        params = self._canonicalize_params(params)
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

    def _flow_key(self, flow_time: float | None, R: int, T: int) -> data_organizer.FlowKey:
        return (data_organizer._normalize_flow_time(flow_time), int(R), int(T))

    def _cache_params_for_wrt(self, flow_time: float | None, R: int, T: int) -> dict[str, Any]:
        params = {"R": int(R), "T": int(T)}
        if flow_time is not None:
            params["flow_time"] = data_organizer._normalize_flow_time(flow_time)
        return params

    def _infer_pair_order(self) -> list[data_organizer.FlowKey]:
        if self._pair_order is not None:
            return self._pair_order

        try:
            obs_L = np.asarray(self.file_data.get("L").values)
            obs_T = np.asarray(self.file_data.get("T").values)
        except ValueError as e:
            raise KeyError(f"Missing required columns (L or T) in file: {e}") from e

        rows_per_cfg = 0
        flow_obs = next((o for o in self.file_data.observables if o.name == "t_over_a2"), None)
        obs_flow = np.asarray(flow_obs.values) if flow_obs is not None else None
        step_obs = next((o for o in self.file_data.observables if o.name in ["# step", "step", "conf_id"]), None)
        if step_obs is not None:
            steps = np.asarray(step_obs.values)
            if len(steps) > 0:
                change_idx = np.flatnonzero(steps != steps[0])
                rows_per_cfg = int(change_idx[0]) if change_idx.size > 0 else int(len(steps))

        if rows_per_cfg <= 0:
            seen: set[data_organizer.FlowKey] = set()
            flow_iter = obs_flow if obs_flow is not None else [None] * len(obs_L)
            for i, (flow_time, l_val, t_val) in enumerate(zip(flow_iter, obs_L, obs_T)):
                pair = self._flow_key(flow_time, int(l_val), int(t_val))
                if pair in seen:
                    rows_per_cfg = i
                    break
                seen.add(pair)
            if rows_per_cfg <= 0:
                rows_per_cfg = int(len(obs_L))

        if obs_flow is None:
            self._pair_order = [
                self._flow_key(None, int(l_val), int(t_val))
                for l_val, t_val in zip(obs_L[:rows_per_cfg], obs_T[:rows_per_cfg])
            ]
        else:
            self._pair_order = [
                self._flow_key(flow_time, int(l_val), int(t_val))
                for flow_time, l_val, t_val in zip(obs_flow[:rows_per_cfg], obs_L[:rows_per_cfg], obs_T[:rows_per_cfg])
            ]
        return self._pair_order

    def _ensure_w_rt_cache(self):
        if self._w_rt_cache is not None:
            return

        try:
            obs_val = np.asarray(self.file_data.get("W_temp").values, dtype=float)
        except ValueError as e:
            raise KeyError(f"Missing required column W_temp in file: {e}") from e

        pair_order = self._infer_pair_order()
        rows_per_cfg = len(pair_order)
        if rows_per_cfg > 0 and len(obs_val) % rows_per_cfg == 0:
            w_matrix = obs_val.reshape(-1, rows_per_cfg)
            self._w_rt_cache = {pair: w_matrix[:, idx] for idx, pair in enumerate(pair_order)}
            return

        try:
            obs_L = np.asarray(self.file_data.get("L").values)
            obs_T = np.asarray(self.file_data.get("T").values)
        except ValueError as e:
            raise KeyError(f"Missing required columns (L or T) in file: {e}") from e

        flow_obs = next((o for o in self.file_data.observables if o.name == "t_over_a2"), None)
        obs_flow = np.asarray(flow_obs.values) if flow_obs is not None else None
        flow_iter = obs_flow if obs_flow is not None else [None] * len(obs_L)
        cache = collections.defaultdict(list)
        for flow_time, l_val, t_val, val in zip(flow_iter, obs_L, obs_T, obs_val):
            cache[self._flow_key(flow_time, int(l_val), int(t_val))].append(val)
        self._w_rt_cache = {k: np.asarray(v, dtype=float) for k, v in cache.items()}
        self._pair_order = sorted(self._w_rt_cache)

    def get_available_flow_times(self) -> list[float | None]:
        self._ensure_w_rt_cache()
        if self._w_rt_cache is None:
            return []
        return sorted(
            {flow_time for flow_time, _, _ in self._w_rt_cache},
            key=lambda value: -1.0 if value is None else float(value),
        )

    def get_available_pairs(self, flow_time: float | None = None) -> list[tuple[int, int]]:
        normalized_flow = data_organizer._normalize_flow_time(flow_time)
        if self._pair_order is not None:
            return [(r, t) for current_flow, r, t in self._pair_order if current_flow == normalized_flow]
        if self._w_rt_cache is not None:
            self._pair_order = sorted(self._w_rt_cache)
            return [(r, t) for current_flow, r, t in self._pair_order if current_flow == normalized_flow]
        return [(r, t) for current_flow, r, t in self._infer_pair_order() if current_flow == normalized_flow]

    def get_unique_Rs(self, r_min: int | None = None, flow_time: float | None = None) -> list[int]:
        rs = sorted({r for r, _ in self.get_available_pairs(flow_time=flow_time)})
        if r_min is not None:
            rs = [r for r in rs if r >= r_min]
        return rs

    def get_unique_Ts(self, R: int | None = None, flow_time: float | None = None) -> list[int]:
        pairs = self.get_available_pairs(flow_time=flow_time)
        if R is None:
            return sorted({t for _, t in pairs})
        return sorted({t for r, t in pairs if r == R})

    def _get_block_layout(self, n_samples: int) -> tuple[np.ndarray, np.ndarray]:
        """
        Partition the ordered Monte Carlo history into contiguous blocks.
        `step_size` is interpreted as the bootstrap block length, not as a
        thinning stride.
        """
        if n_samples <= 0:
            raise ValueError("n_samples must be positive for bootstrap resampling")

        block_starts = np.arange(0, n_samples, self.step_size, dtype=np.int64)
        block_lengths = np.minimum(self.step_size, n_samples - block_starts).astype(np.int64, copy=False)
        return block_starts, block_lengths

    def _get_bootstrap_indices_seq(self, n_samples: int):
        """
        Generate one block-bootstrap replica at a time.
        Each replica resamples contiguous blocks with replacement, so samples
        between block boundaries are reused rather than discarded.
        """
        block_starts, block_lengths = self._get_block_layout(n_samples)
        n_blocks = len(block_starts)
        block_indices = [
            np.arange(start, start + length, dtype=np.int64)
            for start, length in zip(block_starts, block_lengths, strict=False)
        ]
        rng = np.random.default_rng(self.seed)

        for _ in range(self.n_bootstrap):
            sampled_blocks = rng.integers(0, n_blocks, size=n_blocks)
            total_length = int(block_lengths[sampled_blocks].sum())
            indices = np.empty(total_length, dtype=np.int64)
            offset = 0
            for block_idx in sampled_blocks:
                current_block = block_indices[int(block_idx)]
                next_offset = offset + len(current_block)
                indices[offset:next_offset] = current_block
                offset = next_offset
            yield indices

    def _nan_bootstrap_array(self) -> np.ndarray:
        return np.full(self.n_bootstrap, np.nan, dtype=float)

    def _build_compact_block_sums(
        self,
        pairs: list[data_organizer.FlowKey],
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        self._ensure_w_rt_cache()
        if self._w_rt_cache is None:
            raise ValueError("Compact Wilson cache is unavailable.")
        if not pairs:
            return (
                np.array([], dtype=float),
                np.empty((0, 0), dtype=np.float64),
                np.array([], dtype=np.float64),
            )

        first_pair = pairs[0]
        first_values = np.asarray(self._w_rt_cache[first_pair], dtype=float)
        n_samples = int(len(first_values))
        if n_samples <= 0:
            raise ValueError("No Wilson-loop samples available for bootstrap priming.")

        block_starts, block_lengths = self._get_block_layout(n_samples)
        block_sum_matrix = np.empty((len(block_starts), len(pairs)), dtype=np.float64)
        mean_values = np.empty(len(pairs), dtype=float)

        for idx, pair in enumerate(pairs):
            series = self._w_rt_cache.get(pair)
            if series is None:
                raise ValueError(f"No data found for flow_time={pair[0]}, R={pair[1]}, T={pair[2]}")
            series_arr = np.asarray(series, dtype=float)
            if len(series_arr) != n_samples:
                raise ValueError("Compact Wilson data must have a consistent number of samples per pair.")
            mean_values[idx] = float(np.mean(series_arr))
            block_sum_matrix[:, idx] = np.add.reduceat(series_arr.astype(np.float64, copy=False), block_starts)

        return mean_values, block_sum_matrix, block_lengths.astype(np.float64, copy=False)

    def prime_w_rt_cache(
        self,
        pairs: list[tuple[int, int]] | None = None,
        flow_time: float | None = None,
        bootstrap_chunk_size: int = 32,
    ) -> None:
        compact_wilson = getattr(self.file_data, "wilson_by_flow_pair", None)
        if compact_wilson is None:
            return

        if bootstrap_chunk_size < 1:
            raise ValueError("bootstrap_chunk_size must be at least 1")

        normalized_flow = data_organizer._normalize_flow_time(flow_time)
        available_pairs = self.get_available_pairs(flow_time=normalized_flow)
        if pairs is None:
            target_pairs = available_pairs
        else:
            deduped_pairs: list[tuple[int, int]] = []
            seen: set[tuple[int, int]] = set()
            for r_val, t_val in pairs:
                pair = (int(r_val), int(t_val))
                if pair in seen:
                    continue
                deduped_pairs.append(pair)
                seen.add(pair)
            target_pairs = deduped_pairs

        pending_pairs: list[data_organizer.FlowKey] = []
        for pair in target_pairs:
            key = self._flow_key(normalized_flow, pair[0], pair[1])
            if key not in compact_wilson:
                raise ValueError(f"No data found for R={pair[0]}, T={pair[1]}")
            cache_key = make_key("W_R_T", self._cache_params_for_wrt(normalized_flow, pair[0], pair[1]))
            if cache_key not in self.variables:
                pending_pairs.append(key)

        if not pending_pairs:
            return

        pairs_by_length: dict[int, list[data_organizer.FlowKey]] = collections.defaultdict(list)
        for pair in pending_pairs:
            series = compact_wilson[pair]
            pairs_by_length[int(len(series))].append(pair)

        for length, length_pairs in sorted(pairs_by_length.items()):
            if length <= 0:
                continue
            mean_values, block_sum_matrix, block_lengths = self._build_compact_block_sums(length_pairs)
            n_blocks = block_sum_matrix.shape[0]
            bootstrap_matrix = np.empty((self.n_bootstrap, len(length_pairs)), dtype=float)
            rng = np.random.default_rng(self.seed)

            for start in range(0, self.n_bootstrap, bootstrap_chunk_size):
                current_size = min(bootstrap_chunk_size, self.n_bootstrap - start)
                sampled_blocks = rng.integers(0, n_blocks, size=(current_size, n_blocks))
                counts = np.zeros((current_size, n_blocks), dtype=np.float64)
                np.add.at(counts, (np.arange(current_size)[:, None], sampled_blocks), 1.0)
                total_lengths = counts @ block_lengths
                bootstrap_matrix[start:start + current_size] = (
                    (counts @ block_sum_matrix) / total_lengths[:, None]
                ).astype(float, copy=False)

            for idx, pair in enumerate(length_pairs):
                flow_val, r_val, t_val = pair
                cache_key = make_key("W_R_T", self._cache_params_for_wrt(flow_val, r_val, t_val))
                var_data = data_organizer.VariableData("W_R_T")
                var_data.set_value(
                    mean_values[idx],
                    bootstrap_samples=bootstrap_matrix[:, idx],
                    R=r_val,
                    T=t_val,
                    flow_time=flow_val,
                    n_samples=length,
                )
                self.variables[cache_key] = var_data
    
    ### Variable implementations ###

    @register("W_R_T")
    def _calc_W_R_T(self, R: int, T: int, flow_time: float | None = None) -> data_organizer.VariableData: #marker wrt
        key = self._flow_key(flow_time, int(R), int(T))
        self._ensure_w_rt_cache()
        selected_values = self._w_rt_cache.get(key) if self._w_rt_cache is not None else None

        if selected_values is None:
            raise ValueError(f"No data found for flow_time={key[0]}, R={R}, T={T}")

        selected_values = np.asarray(selected_values, dtype=float)
        n_samples = len(selected_values)

        if n_samples == 0:
            raise ValueError(f"No data found for R={R}, T={T}")

        var_data = data_organizer.VariableData("W_R_T")

        mean_val = float(np.mean(selected_values, dtype=np.float64))

        block_starts, block_lengths = self._get_block_layout(n_samples)
        n_blocks = len(block_starts)
        block_sums = np.add.reduceat(
            selected_values.astype(np.float64, copy=False),
            block_starts,
        )

        rng = np.random.default_rng(self.seed)
        sampled_blocks = rng.integers(0, n_blocks, size=(self.n_bootstrap, n_blocks))
        bootstrap_means = (
            block_sums[sampled_blocks].sum(axis=1) /
            block_lengths[sampled_blocks].sum(axis=1)
        ).astype(float, copy=False)

        var_data.set_value(
            mean_val,
            bootstrap_samples=bootstrap_means,
            R=R,
            T=T,
            flow_time=key[0],
            n_samples=n_samples,
        )
        return var_data
    
    @register("V_R")
    def _calc_V_R(
        self,
        R: int,
        t_min: int = 1,
        t_max: int = None,
        flow_time: float | None = None,
    ) -> data_organizer.VariableData:
        """
        Calculates V(R) and its error by fitting W(R,T).
        Includes robust handling for negative signal and array shape mismatches.
        """
        normalized_flow = data_organizer._normalize_flow_time(flow_time)
        unique_ts = self.get_unique_Ts(R=R, flow_time=normalized_flow)

        # Filter T >= t_min and T <= t_max (if provided)
        if t_max is not None:
            ts_to_fit = sorted([t for t in unique_ts if t_min <= t <= t_max])
        else:
            ts_to_fit = sorted([t for t in unique_ts if t >= t_min])
        
        # W_R_T Variables
        w_vars = [
            self.get_variable("W_R_T", **self._cache_params_for_wrt(normalized_flow, R, t))
            for t in ts_to_fit
        ]
        
        ws_main = np.asarray([v.get() for v in w_vars], dtype=float)
        ts_fit  = np.asarray(ts_to_fit, dtype=float)
        ws_err  = np.asarray([v.err() for v in w_vars], dtype=float)

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
            var_data.set_value(np.nan, bootstrap_samples=self._nan_bootstrap_array(), R=R, flow_time=normalized_flow)
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
                return np.nan, np.nan

            try:
                # Log-Linear Fit: ln(W) = ln(C) - V*t
                ln_ws = np.log(ws)
                
                w_weights = None
                if sigma is not None:
                    # Weights for polyfit should be inverse variance of y = ln(W)
                    # var(ln W) ~ (sigma / W)
                    # w = 1/var ~ (W / sigma)
                    # Filter out zero sigma to avoid division by zero (though handled outside)

                    # Attention: np.polyfit squares these weights internally.
                    valid_sigma = (sigma > 0)
                    if np.all(valid_sigma):
                         w_weights = (ws / sigma)
                    else:
                         # Fallback if sigma has zeros (should not happen due to check above)
                         w_weights = None

                # Fit: y = p[0]*x + p[1]  =>  ln(W) = -V*t + ln(C)
                p = np.polyfit(ts, ln_ws, 1, w=w_weights)
                
                V = -p[0]
                C = np.exp(p[1])
                return V, C
            except (RuntimeError, ValueError, TypeError, np.linalg.LinAlgError):
                return np.nan, np.nan

        # Weighted Fit
        V_main, C_main = fit_potential(ts_fit, ws_main, sigma=ws_err)
        
        if np.isnan(V_main):
            # If the main fit fails, return NaN immediately
            var_data = data_organizer.VariableData("V_R")
            var_data.set_value(np.nan, bootstrap_samples=self._nan_bootstrap_array(), R=R, flow_time=normalized_flow)
            return var_data

        # Bootstrap Fits
        n_boot = self.n_bootstrap
        bootstrap_Vs = np.empty(n_boot, dtype=float)

        boot_arrays = [v.bootstrap() for v in w_vars]
        if any(boot is None for boot in boot_arrays):
            raise ValueError(f"Bootstrap samples not available for V_R at R={R}")
        all_boots = np.stack([np.asarray(boot, dtype=float) for boot in boot_arrays], axis=0)
        all_boots_valid = all_boots[:valid_len, :]  # shape: (valid_len, n_boot)
        use_weights = ws_err is not None and np.all(ws_err > 0)

        for i in range(n_boot):
            ws_sample = np.array(all_boots_valid[:, i], copy=True)
            invalid_points = ~np.isfinite(ws_sample) | (ws_sample <= 0)
            if np.any(invalid_points):
                ws_sample[invalid_points] = ws_main[invalid_points]

            if np.any(ws_sample <= 0):
                bootstrap_Vs[i] = V_main
                continue

            ln_ws_sample = np.log(ws_sample)
            w_weights = (ws_sample / ws_err) if use_weights else None

            try:
                p = np.polyfit(ts_fit, ln_ws_sample, 1, w=w_weights)
                bootstrap_Vs[i] = -p[0]
            except (RuntimeError, ValueError, TypeError, np.linalg.LinAlgError):
                bootstrap_Vs[i] = V_main

        var_data = data_organizer.VariableData("V_R")
        var_data.set_value(
            V_main,
            bootstrap_samples=bootstrap_Vs,
            R=R,
            t_min=t_min,
            t_max=t_max,
            fit_C=C_main,
            flow_time=normalized_flow,
        )
        return var_data
    
    @register("tau_int")
    def _calc_tau_int(self, obs_name: str = "plaquette", S: float = 1.5) -> data_organizer.VariableData:
        if obs_name == "W_temp":
            raise KeyError("Cannot run FFT Autocorrelation on interleaved W_temp data.")

        try:
            obs = self.get_observable(obs_name)
        except KeyError:
            # Fallback for commonly used names if exact match fails
            found = False
            for alias in ["plaquette", "retrace"]:
                 for o in self.file_data.observables:
                     if o.name == alias:
                         obs = o
                         found = True
                         break
                 if found:
                     break
            if not found:
                 var_data = data_organizer.VariableData("tau_int")
                 var_data.set_value(0.5, obs_name=obs_name)
                 return var_data

        series = np.asarray(obs.values)
        n = len(series)
        
        # Default fallback if series is too short
        tau = 0.5 

        if n >= 100:
            # Subtract mean
            series_centered = series - np.mean(series)
            var = np.var(series_centered)
            
            if var > 0:
                # FFT for autocorrelation
                ft = np.fft.rfft(series_centered, n=2*n)
                spec = np.abs(ft)**2
                acf = np.fft.irfft(spec, n=2*n)
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
    def _calc_r0(
        self,
        t_min: int = 5,
        t_max: int = None,
        target_force: float = 1.65,
        r_min: int = 2,
        flow_time: float | None = None,
    ) -> data_organizer.VariableData:
        """
        Calculates the Sommer parameter r0/a by fitting the Cornell potential to V(R).
        Replaces analyze_wilson.fit_sommer_parameter.
        """
        def format_point_diagnostics(point_rows):
            if not point_rows:
                return "no V(R) candidates were collected"
            return "; ".join(
                f"R={row['R']}: V={row['V']}, err={row['err']}, boot_finite={row['boot_finite']}/{row['boot_total']} -> {row['status']}"
                for row in point_rows
            )

        # 1. Identify available R values
        normalized_flow = data_organizer._normalize_flow_time(flow_time)
        try:
            unique_Rs = self.get_unique_Rs(r_min=r_min, flow_time=normalized_flow)
        except KeyError:
             raise KeyError("Could not determine available R (L) values from file data.")

        # 2. Gather V(R) for all R
        rs = []
        vs = []
        errs = []
        point_diagnostics = []
        
        # We need a consistent set of bootstrap samples to propagate errors correctly
        # So we collect the bootstrap arrays from each V_R
        bootstrap_rows: list[np.ndarray] = []
        for r in unique_Rs:
            try:
                # Reuse the existing V_R calculator
                v_var = self.get_variable("V_R", R=int(r), t_min=t_min, t_max=t_max, flow_time=normalized_flow)
                
                val = v_var.get()
                boot_samples = v_var.bootstrap()
                boot_arr = np.asarray(boot_samples, dtype=float) if boot_samples is not None else np.array([], dtype=float)
                finite_boot_count = int(np.isfinite(boot_arr).sum())
                err_val = v_var.err()
                if np.isnan(val) or np.isinf(val):
                    point_diagnostics.append({
                        "R": int(r),
                        "V": val,
                        "err": err_val,
                        "boot_finite": finite_boot_count,
                        "boot_total": int(boot_arr.size),
                        "status": "discarded (non-finite V)",
                    })
                    continue

                rs.append(r)
                vs.append(val)
                errs.append(err_val if err_val is not None else 1.0)
                point_diagnostics.append({
                    "R": int(r),
                    "V": val,
                    "err": err_val,
                    "boot_finite": finite_boot_count,
                    "boot_total": int(boot_arr.size),
                    "status": "candidate",
                })
                
                # Collect bootstraps. Must ensure they align (same seed/length).
                # The Calculator class ensures consistent seeding in __init__.
                if boot_samples is not None:
                    bootstrap_rows.append(boot_arr)
                
            except (ValueError, RuntimeError) as exc:
                point_diagnostics.append({
                    "R": int(r),
                    "V": np.nan,
                    "err": None,
                    "boot_finite": 0,
                    "boot_total": 0,
                    "status": f"discarded ({type(exc).__name__}: {exc})",
                })
                continue
        if len(rs) < 3:
            raise ValueError(
                f"Not enough valid V(R) points to fit r0 (found {len(rs)}). "
                f"Diagnostics: {format_point_diagnostics(point_diagnostics)}"
            )

        rs_arr = np.asarray(rs, dtype=float)
        vs_arr = np.asarray(vs, dtype=float)
        errs_arr = np.asarray(errs, dtype=float)
        boot_matrix = np.stack(bootstrap_rows, axis=0) if len(bootstrap_rows) == len(rs_arr) and bootstrap_rows else None

        try:
            fit_result = fit_r0_from_potential_data(
                rs_arr,
                vs_arr,
                errs=errs_arr,
                bootstrap_matrix=boot_matrix,
                target_force=target_force,
                n_bootstrap=self.n_bootstrap,
            )
        except ValueError as exc:
            raise ValueError(
                f"{exc} Diagnostics: {format_point_diagnostics(point_diagnostics)}"
            ) from exc

        var_data = data_organizer.VariableData("r0")
        var_data.set_value(
            fit_result["r0"],
            bootstrap_samples=fit_result["bootstrap_samples"],
            t_min=t_min,
            t_max=t_max,
            r_min=r_min,
            cornell_params=fit_result["cornell_params"],
            flow_time=normalized_flow,
        )
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
    def _calc_chi(
        self,
        R: float,
        T: float,
        flow_time: float | None = None,
    ) -> data_organizer.VariableData:
        """
        Calculates the Creutz ratio chi(R, T) using Wilson loops.
        R and T are the half-integer coordinates of the center of the Creutz plaquette.
        e.g., R=1.5, T=1.5 for the ratio of W(1,1), W(2,1), W(1,2), W(2,2).
        """
        # Creutz ratio is defined at half-integer coordinates
        # chi(t + a/2, r + a/2), so r_int = r - 0.5, t_int = t - 0.5
        r_int = int(R - 0.5)
        t_int = int(T - 0.5)
        normalized_flow = data_organizer._normalize_flow_time(flow_time)
        
        try:
            # Get the four Wilson loops needed for the ratio
            w_tr = self.get_variable("W_R_T", **self._cache_params_for_wrt(normalized_flow, r_int, t_int))
            w_t1_r = self.get_variable("W_R_T", **self._cache_params_for_wrt(normalized_flow, r_int, t_int + 1))
            w_t_r1 = self.get_variable("W_R_T", **self._cache_params_for_wrt(normalized_flow, r_int + 1, t_int))
            w_t1_r1 = self.get_variable("W_R_T", **self._cache_params_for_wrt(normalized_flow, r_int + 1, t_int + 1))
        except (ValueError, KeyError) as e:
            raise ValueError(f"Could not get required Wilson loops for chi({R}, {T}): {e}")

        # Main value
        chi_val = creutz_chi_from_wilson(
            w_tr.get(),
            w_t1_r.get(),
            w_t_r1.get(),
            w_t1_r1.get(),
        )
        
        # Bootstrap values
        # Get bootstrap samples for all 4 loops
        b_t1_r = w_t1_r.bootstrap()
        b_t_r1 = w_t_r1.bootstrap()
        b_tr   = w_tr.bootstrap()
        b_t1_r1= w_t1_r1.bootstrap()

        # Check if bootstrap samples are available
        if b_t1_r is None or b_t_r1 is None or b_tr is None or b_t1_r1 is None:
             raise ValueError(f"Bootstrap samples not available for one of the Wilson loops for chi({R},{T})")

        chi_boots = creutz_chi_from_wilson(b_tr, b_t1_r, b_t_r1, b_t1_r1)
        
        var_data = data_organizer.VariableData("chi")
        var_data.set_value(chi_val, bootstrap_samples=chi_boots, R=R, T=T, flow_time=normalized_flow)
        return var_data

    @register("F_chi")
    def _calc_F_chi(
        self,
        R: float,
        t_large: int,
        flow_time: float | None = None,
    ) -> data_organizer.VariableData:
        """
        Extracts the static force F(R) from Creutz ratios chi(R, T)
        by taking the value at a large time extent, T = t_large.
        R is a half-integer (e.g., 1.5, 2.5).
        t_large is an integer, so T for chi will be t_large + 0.5.
        """
        T = float(t_large) + 0.5
        normalized_flow = data_organizer._normalize_flow_time(flow_time)
        try:
            chi_var = self.get_variable("chi", R=R, T=T, flow_time=normalized_flow)
            var_data = data_organizer.VariableData("F_chi")
            # The value and bootstrap samples are directly from chi(R, t_large)
            var_data.set_value(
                chi_var.get(),
                bootstrap_samples=chi_var.bootstrap(),
                R=R,
                t_large=t_large,
                flow_time=normalized_flow,
            )
            return var_data
        except (ValueError, KeyError) as e:
            raise ValueError(f"Could not calculate F_chi for R={R} at t_large={t_large}: {e}")

    @register("r0_chi")
    def _calc_r0_chi(self, t_large: int = 4, target_force: float = 1.65, max_rel_err: float = 0.5, use_weighted_fit: bool = True, fit_window: int = 2, discard_negative: bool = True, r_min: int = 1, flow_time: float | None = None) -> data_organizer.VariableData:
        """
        Calculates Sommer parameter r0/a from Creutz ratios.
        It interpolates r^2 * F(r) to find where it equals target_force.
        Points with a high relative error on F(r) can be filtered out.
        A weighted linear fit over a window of points can be used for interpolation.
        
        Args:
            t_large: The time extent to use for F_chi extraction.
            target_force: The target value for r^2 * F(r) (default 1.65).
            max_rel_err: Filter out F_chi points with relative error > max_rel_err.
            use_weighted_fit: Use error-weighted fitting for finding r0.
            fit_window: Window size for local fitting.
            discard_negative: If True, discard negative F_chi values before solving.
            r_min: Minimum R value to include in the fit (default 1).
        """
        # 1. Identify available R values from data.
        normalized_flow = data_organizer._normalize_flow_time(flow_time)
        try:
            unique_Ls = self.get_unique_Rs(r_min=r_min, flow_time=normalized_flow)
            unique_Rs = [float(L) + 0.5 for L in unique_Ls if L + 1 in unique_Ls]
        except KeyError:
             raise KeyError("Could not determine available R (L) values from file data for r0_chi.")

        if not unique_Rs:
            raise ValueError("Not enough adjacent L values for standard Creutz ratios.")

        # 2. Gather F_chi(R) for all available R at a fixed large T.
        rs, fs, errs, bootstrap_rows = [], [], [], []
        
        for r_val in unique_Rs:
            try:
                f_var = self.get_variable("F_chi", R=r_val, t_large=t_large, flow_time=normalized_flow)
                val, err = f_var.get(), f_var.err()

                if np.isnan(val) or np.isinf(val) or np.isnan(err) or np.isinf(err):
                    continue
                
                # Filter negative values if requested
                if discard_negative and val <= 0:
                    continue

                # Filter points with large relative error
                if max_rel_err is not None and val != 0 and abs(err / val) > max_rel_err:
                    continue

                rs.append(r_val)
                fs.append(val)
                errs.append(err)
                if f_var.bootstrap() is not None:
                    bootstrap_rows.append(np.asarray(f_var.bootstrap(), dtype=float))
            except (ValueError, KeyError):
                continue

        if len(rs) < 2: # Need at least 2 points for interpolation/fit
            # If we don't have enough points, we can't do anything.
            # But returning a VariableData with None value is safer than crashing
            # if we just filtered everything out.
            var_data = data_organizer.VariableData("r0_chi")
            var_data.set_value(None, t_large=t_large, flow_time=normalized_flow)
            return var_data
        
        rs = np.asarray(rs, dtype=float)
        fs = np.asarray(fs, dtype=float)
        errs = np.asarray(errs, dtype=float)

        def get_r0_from_force(r_vals, f_vals, f_errs=None):
            """
            Finds r0 by linearly interpolating/fitting r^2*F(r).
            Uses global spline interpolation if possible, fallback to local fits.
            """
            if len(r_vals) < 2: return np.nan
            
            r2_vals = r_vals**2
            r2f_vals = r2_vals * f_vals

            # A. Global spline interpolation
            if use_weighted_fit and f_errs is not None and len(f_errs) == len(r_vals):
                from scipy.interpolate import UnivariateSpline
                std_devs = r2_vals * f_errs
                valid = std_devs > 0
                if np.sum(valid) >= 4:
                    r2_v = r2_vals[valid]
                    y_v = r2f_vals[valid] - target_force
                    w_v = 1.0 / std_devs[valid]
                    
                    # Ensure x is strictly increasing
                    sort_idx = np.argsort(r2_v)
                    r2_v = r2_v[sort_idx]
                    y_v = y_v[sort_idx]
                    w_v = w_v[sort_idx]
                    
                    try:
                        spline = UnivariateSpline(r2_v, y_v, w=w_v, k=3)
                        roots = spline.roots()
                        positive_roots = [r for r in roots if r > 0]
                        if positive_roots:
                            return np.sqrt(positive_roots[0])
                    except Exception:
                        pass # Fallback to local fit

            # B. Local crossover fallback
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

                        # Attention: np.polyfit squares these weights internally.
                        win_r4f_err = (win_r2) * (win_f_errs)
                        valid_w = win_r4f_err > 0
                        if np.sum(valid_w) < 2: continue # Not enough points with valid weights

                        weights = np.zeros_like(win_r4f_err)
                        weights[valid_w] = 1.0 / win_r4f_err[valid_w]

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
            
            # B. If no crossover found, and we want to fit globally (or robustly)
            if use_weighted_fit and f_errs is not None and len(f_errs) == len(r_vals):
                 # Fit a line to r^2*F(r) vs r^2 for all points?
                 # Or maybe r^2*F(r) vs r?
                 # Assuming approximate linearity of r^2 F(r) vs r^2 (string tension dominance)
                 # F(r) ~ sigma => r^2 F(r) ~ sigma * r^2.
                 
                 x_fit = r2_vals
                 y_fit = r2f_vals
                 
                 weights_fit = np.zeros_like(x_fit)
                 var_y = (x_fit * f_errs)**2
                 valid_w = var_y > 0
                 if np.sum(valid_w) >= 2:
                     weights_fit[valid_w] = 1.0 / var_y[valid_w]
                     try:
                         # Linear fit: y = slope * x + intercept
                         # r^2 F = slope * r^2 + intercept
                         p = np.polyfit(x_fit[valid_w], y_fit[valid_w], 1, w=weights_fit[valid_w])
                         slope, intercept = p[0], p[1]
                         
                         # Solve slope * r0^2 + intercept = target_force
                         if slope != 0:
                            r0_sq = (target_force - intercept) / slope
                            if r0_sq > 0: return np.sqrt(r0_sq)
                     except (np.linalg.LinAlgError, ValueError):
                         pass

            return np.nan # No solution found

        # Main value
        r0_val = get_r0_from_force(rs, fs, errs)

        # Bootstrap
        fill_value = float(r0_val) if r0_val is not None and np.isfinite(r0_val) else np.nan
        r0_bootstraps = self._nan_bootstrap_array()
        if len(bootstrap_rows) > 0 and len(bootstrap_rows) == len(rs):
            boot_data = np.stack(bootstrap_rows, axis=0).T
            r0_bootstraps.fill(fill_value)
            for i, sample_fs in enumerate(boot_data):
                repaired_fs = np.array(sample_fs, copy=True)
                invalid_points = ~np.isfinite(repaired_fs)
                if np.any(invalid_points):
                    repaired_fs[invalid_points] = fs[invalid_points]

                # For bootstrap samples, we don't have individual errors,
                # so we use the errors from the main data set as weights.
                val = get_r0_from_force(rs, repaired_fs, errs)
                r0_bootstraps[i] = fill_value if np.isnan(val) else val
        elif r0_val is not None and np.isfinite(r0_val):
            r0_bootstraps = np.full(self.n_bootstrap, r0_val, dtype=float)

        var_data = data_organizer.VariableData("r0_chi")
        var_data.set_value(r0_val, bootstrap_samples=r0_bootstraps, t_large=t_large, r_min=r_min, flow_time=normalized_flow)
        return var_data
    
    @register("creutz_P")
    def _calc_creutz_P(self, R: int, flow_time: float | None = None) -> data_organizer.VariableData:
        """
        Calculates the physical ratio P(R) defined by Creutz (1981):
        P(R) = 1 - [W(R,R) * W(R/2,R/2)] / [W(R,R/2)]^2
        
        This quantity is used to study RG flow and remove UV divergences.
        R must be an even integer (representing the larger scale r/a).
        """
        if R % 2 != 0:
            raise ValueError(f"R must be even for Creutz P-ratio calculation (got {R}).")
        
        R_half = R // 2
        normalized_flow = data_organizer._normalize_flow_time(flow_time)
        
        try:
            # Retrieve required Wilson loops: W(R,R), W(R/2,R/2), and W(R, R/2)
            w_large = self.get_variable("W_R_T", **self._cache_params_for_wrt(normalized_flow, R, R))
            w_small = self.get_variable("W_R_T", **self._cache_params_for_wrt(normalized_flow, R_half, R_half))
            
            # Try getting rectangular loop W(R, R/2). 
            # If not found, check for the symmetric W(R/2, R).
            try:
                w_rect = self.get_variable("W_R_T", **self._cache_params_for_wrt(normalized_flow, R, R_half))
            except (ValueError, KeyError):
                w_rect = self.get_variable("W_R_T", **self._cache_params_for_wrt(normalized_flow, R_half, R))
                
        except (ValueError, KeyError) as e:
            raise ValueError(f"Required Wilson loops for creutz_P({R}) not found: {e}")

        def calculate_P(w_LL, w_ss, w_Ls):
            # w_LL = W(R,R), w_ss = W(R/2,R/2), w_Ls = W(R,R/2)
            if w_Ls == 0:
                return np.nan
            
            numerator = w_LL * w_ss
            denominator = w_Ls ** 2
            
            # P = 1 - (W(R,R)W(R/2,R/2)) / W(R,R/2)^2
            return 1.0 - (numerator / denominator)  # Attention: Is 1- necesssary?

        def calculate_P_array(w_LL, w_ss, w_Ls):
            w_ll_arr = np.asarray(w_LL, dtype=float)
            w_ss_arr = np.asarray(w_ss, dtype=float)
            w_ls_arr = np.asarray(w_Ls, dtype=float)
            numerator = w_ll_arr * w_ss_arr
            denominator = w_ls_arr ** 2
            result = np.full(numerator.shape, np.nan, dtype=float)
            valid = np.isfinite(numerator) & np.isfinite(denominator) & (w_ls_arr != 0)
            result[valid] = 1.0 - (numerator[valid] / denominator[valid])
            return result

        # 1. Calculate Main Value
        val_P = calculate_P(w_large.get(), w_small.get(), w_rect.get())
        
        # 2. Calculate Bootstrap Values
        b_large = w_large.bootstrap()
        b_small = w_small.bootstrap()
        b_rect  = w_rect.bootstrap()
        
        # Ensure we have bootstraps for all components
        if b_large is not None and b_small is not None and b_rect is not None:
            P_boots = calculate_P_array(b_large, b_small, b_rect)
        else:
            # If bootstraps are missing (e.g. from a fallback), fill with NaNs or empty
            P_boots = self._nan_bootstrap_array()

        var_data = data_organizer.VariableData("creutz_P")
        var_data.set_value(val_P, bootstrap_samples=P_boots, R=R, flow_time=normalized_flow)
        return var_data
    
    @register("a_creutz")
    def _calc_a_creutz(self, R: int, sigma_sqrt_MeV: float = 440.0, flow_time: float | None = None) -> data_organizer.VariableData:
        """
        Calculates lattice spacing 'a' in fm from Creutz ratio P(R).
        Assumes Area Law dominance: P(R) = 1 - exp(-1/4 * sigma * a^2 * R^2).
        sigma_sqrt_MeV: Physical string tension sqrt(sigma) in MeV (default 440).
        """
        normalized_flow = data_organizer._normalize_flow_time(flow_time)
        try:
            p_var = self.get_variable("creutz_P", R=R, flow_time=normalized_flow)
        except (ValueError, KeyError) as e:
             raise ValueError(f"Could not calculate a_creutz({R}): {e}")

        hc = 197.3 # MeV fm

        def calculate_a(p_val):
            # p_val must be < 1.0 for real result. 
            # In confinement region, p_val approaches 1.0 from below.
            # Also if p_val <= 0, log is undefined (or imaginary).
            if np.isnan(p_val) or p_val >= 1.0 or p_val <= 0.0: 
                return np.nan
            
            term = np.sqrt(-np.log(1.0 - p_val))
            # a = (2 / (R * sqrt(sigma))) * term
            # Convert to fm using hc / sigma_sqrt
            a_fm = (hc / sigma_sqrt_MeV) * (2.0 / R) * term
            return a_fm

        def calculate_a_array(p_vals):
            p_arr = np.asarray(p_vals, dtype=float)
            result = np.full(p_arr.shape, np.nan, dtype=float)
            valid = np.isfinite(p_arr) & (p_arr < 1.0) & (p_arr > 0.0)
            result[valid] = (hc / sigma_sqrt_MeV) * (2.0 / R) * np.sqrt(-np.log(1.0 - p_arr[valid]))
            return result

        # Main value
        val_a = calculate_a(p_var.get())
        
        # Bootstrap
        p_boots = p_var.bootstrap()
        
        if p_boots is not None:
            a_boots = calculate_a_array(p_boots)
        else:
             a_boots = self._nan_bootstrap_array()
             
        var_data = data_organizer.VariableData("a_creutz")
        var_data.set_value(val_a, bootstrap_samples=a_boots, R=R, sigma_sqrt=sigma_sqrt_MeV, flow_time=normalized_flow)
        return var_data
    
    @register("volume_r0")
    def _calc_volume(
        self,
        L: int | None = None,
        L0: int | None = None,
        L1: int | None = None,
        L2: int | None = None,
        L3: int | None = None,
        r0_phys: float = 0.5,
        **r0_kwargs,
    ) -> data_organizer.VariableData:
        # Physical 4D volume based on r0. For isotropic legacy calls this reduces to
        # (L*a)^4; otherwise use (L0*L1*L2*L3)*a^4.
        
        try:
            r0_var = self.get_variable("r0", **r0_kwargs)
        except (ValueError, KeyError) as e:
             raise ValueError(f"Could not calculate r0 for volume: {e}")

        r0_on_a = r0_var.get()

        if all(val is not None for val in (L0, L1, L2, L3)):
            lattice_extent = int(L0) * int(L1) * int(L2) * int(L3)
        elif L is not None:
            lattice_extent = int(L) ** 4
        else:
            raise ValueError("volume_r0 requires either L or all of L0, L1, L2, and L3")

        def calculate_volume(val_r0_on_a):
            if val_r0_on_a is None or np.isnan(val_r0_on_a) or val_r0_on_a <= 0:
                return np.nan
            return calculate_volume_from_r0(lattice_extent, val_r0_on_a, r0_phys)

        def calculate_volume_array(values):
            arr = np.asarray(values, dtype=float)
            result = np.full(arr.shape, np.nan, dtype=float)
            valid = np.isfinite(arr) & (arr > 0)
            result[valid] = lattice_extent * ((r0_phys / arr[valid]) ** 4)
            return result

        # Main Value
        vol_val = calculate_volume(r0_on_a)
        
        # Bootstrap
        r0_boots = r0_var.bootstrap()
        
        if r0_boots is not None:
             vol_boots = calculate_volume_array(r0_boots)
        else:
             vol_boots = self._nan_bootstrap_array()
             
        var_data = data_organizer.VariableData("volume_r0")
        var_data.set_value(
            vol_val,
            bootstrap_samples=vol_boots,
            L=L,
            L0=L0,
            L1=L1,
            L2=L2,
            L3=L3,
            r0_phys=r0_phys,
        )
        return var_data

    @register("epsilon_bar")
    def _calc_epsilon_bar(
        self,
        epsilon1: float,
        beta: float,
        **r0_kwargs,
    ) -> data_organizer.VariableData:
        if epsilon1 == 0:
            zero_boots = np.zeros(self.n_bootstrap, dtype=float)
            var_data = data_organizer.VariableData("epsilon_bar")
            var_data.set_value(
                0.0,
                bootstrap_samples=zero_boots,
                epsilon1=epsilon1,
                beta=beta,
                **r0_kwargs,
            )
            return var_data

        try:
            r0_var = self.get_variable("r0", **r0_kwargs)
        except (ValueError, KeyError) as e:
             raise ValueError(f"Could not calculate r0 for epsilon_bar: {e}")

        def calculate_epsilon_bar(val_r0_on_a):
            if beta == 0 or val_r0_on_a is None or np.isnan(val_r0_on_a):
                return np.nan
            # Here r0 denotes the fitted Sommer scale r0/a, so
            # epsilon_bar = epsilon1 / beta * (r0 / a)^2 = epsilon1 / beta * (r0/a)^2.
            return (epsilon1 / beta) * (val_r0_on_a ** 2)

        def calculate_epsilon_bar_array(values):
            arr = np.asarray(values, dtype=float)
            result = np.full(arr.shape, np.nan, dtype=float)
            if beta == 0:
                return result
            valid = np.isfinite(arr)
            result[valid] = (epsilon1 / beta) * (arr[valid] ** 2)
            return result

        eps_bar_val = calculate_epsilon_bar(r0_var.get())

        r0_boots = r0_var.bootstrap()
        if r0_boots is not None:
            eps_bar_boots = calculate_epsilon_bar_array(r0_boots)
        else:
            eps_bar_boots = self._nan_bootstrap_array()

        var_data = data_organizer.VariableData("epsilon_bar")
        var_data.set_value(
            eps_bar_val,
            bootstrap_samples=eps_bar_boots,
            epsilon1=epsilon1,
            beta=beta,
            **r0_kwargs,
        )
        return var_data

    @register("eps_bar")
    def _calc_eps_bar(
        self,
        epsilon1: float,
        beta: float,
        **r0_kwargs,
    ) -> data_organizer.VariableData:
        return self._calc_epsilon_bar(epsilon1=epsilon1, beta=beta, **r0_kwargs)

    def _calculate_volume_value(self, lattice_extent, val_r0_on_a, r0_phys):
        return calculate_volume_from_r0(lattice_extent, val_r0_on_a, r0_phys)
    
    @register("length")
    def _calc_length(
        self,
        L0: int,
        r0_phys: float = 0.5,
        **r0_kwargs,
    ) -> data_organizer.VariableData:
        try:
            r0_var = self.get_variable("r0", **r0_kwargs)
        except (ValueError, KeyError) as e:
             raise ValueError(f"Could not calculate r0 for length: {e}")

        r0_on_a = r0_var.get()

        def calculate_length(val_r0_on_a):
            if val_r0_on_a is None or np.isnan(val_r0_on_a) or val_r0_on_a <= 0:
                return np.nan
            return L0 * (r0_phys / val_r0_on_a)

        def calculate_length_array(values):
            arr = np.asarray(values, dtype=float)
            result = np.full(arr.shape, np.nan, dtype=float)
            valid = np.isfinite(arr) & (arr > 0)
            result[valid] = L0 * (r0_phys / arr[valid])
            return result

        # Main Value
        len_val = calculate_length(r0_on_a)
        
        # Bootstrap
        r0_boots = r0_var.bootstrap()
        
        if r0_boots is not None:
             len_boots = calculate_length_array(r0_boots)
        else:
             len_boots = self._nan_bootstrap_array()
             
        var_data = data_organizer.VariableData("length")
        var_data.set_value(
            len_val,
            bootstrap_samples=len_boots,
            L0=L0,
            r0_phys=r0_phys,
        )
        return var_data

    @register("effective_mass")
    def _calc_effective_mass(
        self,
        R: int,
        T: int = 1,
        flow_time: float | None = None,
    ) -> data_organizer.VariableData:
        """
        Calculates the effective mass m_eff(R) from Wilson loops W(R,T).
        m_eff(R) = ln[W(R,T)/W(R,T+1)] for a fixed R and large T.
        """
        normalized_flow = data_organizer._normalize_flow_time(flow_time)
        try:
            # Get W(R,T) for T and T+1
            w_t = self.get_variable("W_R_T", **self._cache_params_for_wrt(normalized_flow, R, T))
            w_t1 = self.get_variable("W_R_T", **self._cache_params_for_wrt(normalized_flow, R, T + 1))
        except (ValueError, KeyError) as e:
            raise ValueError(f"Could not get required Wilson loops for effective mass at R={R}: {e}")

        def calculate_m_eff(w_T, w_T1):
            if w_T1 == 0 or w_T <= 0 or w_T1 <= 0:
                return np.nan
            return np.log(w_T / w_T1)

        def calculate_m_eff_array(w_T, w_T1):
            w_t_arr = np.asarray(w_T, dtype=float)
            w_t1_arr = np.asarray(w_T1, dtype=float)
            result = np.full(w_t_arr.shape, np.nan, dtype=float)
            valid = np.isfinite(w_t_arr) & np.isfinite(w_t1_arr) & (w_t_arr > 0) & (w_t1_arr > 0)
            result[valid] = np.log(w_t_arr[valid] / w_t1_arr[valid])
            return result

        # Main value
        m_eff_val = calculate_m_eff(w_t.get(), w_t1.get())
        
        # Bootstrap values
        b_t = w_t.bootstrap()
        b_t1 = w_t1.bootstrap()

        if b_t is not None and b_t1 is not None:
            m_eff_boots = calculate_m_eff_array(b_t, b_t1)
        else:
             m_eff_boots = self._nan_bootstrap_array()

        var_data = data_organizer.VariableData("effective_mass")
        var_data.set_value(m_eff_val, bootstrap_samples=m_eff_boots, R=R, T=T, flow_time=normalized_flow)
        return var_data
