import json
import hashlib
import os
import threading
import concurrent.futures
from dataclasses import asdict
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple

import numpy as np

from load_input_yaml import (
    GaugeObservableParams,
    GradientFlowParams,
    MetropolisParams,
    load_gradient_flow_params,
    load_params,
)

import data_organizer as do
from calculator import Calculator, make_key

# --- CONFIGURATION ---
DATA_ROOT = Path("../data").resolve()
CALC_RESULT_BASE = DATA_ROOT / "calcResult"
THERMALIZATION_STEPS = 1500
CALC_VERSION = "5.1"
GROUP_IGNORE_METRO_FIELDS = {"seed", "nSweep"}

DEFAULT_N_BOOTSTRAP = 200
DEFAULT_BLOCK_SIZE = 2500
DEFAULT_SOMMER_TARGET = 1.65
DEFAULT_V_R_T_MIN = 1
DEFAULT_V_R_T_MAX = None
DEFAULT_R0_T_MIN = 8
DEFAULT_R0_T_MAX = None
DEFAULT_R0_R_MIN = 2
DEFAULT_R0_CHI_T_LARGE = 4
DEFAULT_R0_CHI_MAX_REL_ERR = 0.5
DEFAULT_R0_CHI_USE_WEIGHTED_FIT = True
DEFAULT_R0_CHI_FIT_WINDOW = 2
DEFAULT_R0_CHI_DISCARD_NEGATIVE = True
DEFAULT_R0_CHI_R_MIN = 1

_GROUP_INDEX_CACHE: Dict[str, List[str]] | None = None
_GROUP_INDEX_ROOT: Path | None = None


def _default_analysis_options() -> Dict[str, Any]:
    return {
        "r0_t_min": DEFAULT_R0_T_MIN,
    }


def _resolve_analysis_options(analysis_options: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    resolved = _default_analysis_options()
    if analysis_options:
        resolved.update(analysis_options)

    resolved["r0_t_min"] = int(resolved["r0_t_min"])
    return resolved


def _cache_id_for_analysis(base_id: str, analysis_options: Dict[str, Any]) -> str:
    defaults = _default_analysis_options()
    if analysis_options == defaults:
        return base_id

    payload = json.dumps(analysis_options, sort_keys=True, separators=(",", ":"))
    digest = hashlib.md5(payload.encode("utf-8")).hexdigest()[:12]
    return f"{base_id}__cfg_{digest}"


def _range_summary(values: List[float | int]) -> Dict[str, Any]:
    if not values:
        return {"min": None, "max": None, "count": 0}
    return {
        "min": float(values[0]) if isinstance(values[0], float) else int(values[0]),
        "max": float(values[-1]) if isinstance(values[-1], float) else int(values[-1]),
        "count": len(values),
    }


def _build_analysis_settings(
    unique_Rs: List[int],
    unique_Ts: List[int],
    tau: float,
    block_size: int,
    sommer_target: float,
    r0_t_min: int,
) -> Dict[str, Any]:
    unique_R_set = set(unique_Rs)
    r0_fit_Rs = [int(r) for r in unique_Rs if int(r) >= DEFAULT_R0_R_MIN]
    r0_chi_fit_Rs = [r + 0.5 for r in unique_Rs if (r + 1) in unique_R_set]

    return {
        "thermalization_steps": THERMALIZATION_STEPS,
        "sommer_target": float(sommer_target),
        "n_bootstrap": DEFAULT_N_BOOTSTRAP,
        "tau_int": float(tau),
        "bootstrap_block_size": int(block_size),
        "available_R": _range_summary(unique_Rs),
        "available_T": _range_summary(unique_Ts),
        "V_R": {
            "t_min": DEFAULT_V_R_T_MIN,
            "t_max": DEFAULT_V_R_T_MAX,
            "r_min": min(unique_Rs) if unique_Rs else None,
            "r_max": max(unique_Rs) if unique_Rs else None,
        },
        "r0": {
            "t_min": int(r0_t_min),
            "t_max": DEFAULT_R0_T_MAX,
            "r_min": DEFAULT_R0_R_MIN,
            "r_max": max(r0_fit_Rs) if r0_fit_Rs else None,
        },
        "r0_chi": {
            "t_large": DEFAULT_R0_CHI_T_LARGE,
            "r_min": min(r0_chi_fit_Rs) if r0_chi_fit_Rs else None,
            "r_max": max(r0_chi_fit_Rs) if r0_chi_fit_Rs else None,
            "max_rel_err": DEFAULT_R0_CHI_MAX_REL_ERR,
            "use_weighted_fit": DEFAULT_R0_CHI_USE_WEIGHTED_FIT,
            "fit_window": DEFAULT_R0_CHI_FIT_WINDOW,
            "discard_negative": DEFAULT_R0_CHI_DISCARD_NEGATIVE,
        },
        "creutz_ratio": {
            "definition": "log(W(R,T+1)*W(R+1,T)/(W(R,T)*W(R+1,T+1)))",
            "estimator": "ratio of Wilson-loop means with block-bootstrap errors",
        },
    }


def get_run_id(path: str) -> str:
    path_abs = os.path.abspath(path)
    rel_path = os.path.relpath(path_abs, start=str(DATA_ROOT))
    return hashlib.md5(rel_path.encode("utf-8")).hexdigest()


def get_result_path(cache_id: str) -> Path:
    return CALC_RESULT_BASE / f"{cache_id}.json"


def load_cached_result(cache_id: str) -> Optional[Dict[str, Any]]:
    p = get_result_path(cache_id)
    if p.exists():
        try:
            with open(p, "r") as f:
                data = json.load(f)
                if data.get("version") == CALC_VERSION:
                    return data
        except (json.JSONDecodeError, IOError):
            return None
    return None


def save_result(cache_id: str, data: Dict[str, Any]):
    CALC_RESULT_BASE.mkdir(parents=True, exist_ok=True)
    data["version"] = CALC_VERSION
    p = get_result_path(cache_id)
    with open(p, "w") as f:
        json.dump(data, f, indent=2)


def _group_key_from_params(
    metro: MetropolisParams,
    gauge: GaugeObservableParams,
    gradient_flow: GradientFlowParams | None = None,
) -> str:
    metro_dict = asdict(metro)
    for key in GROUP_IGNORE_METRO_FIELDS:
        metro_dict.pop(key, None)
    if gradient_flow is None:
        gradient_flow = GradientFlowParams(
            enabled=False,
            integrator="",
            dt=0.0,
            t_values=[],
            measure_energy_clover=False,
            measure_wilson_loop_temporal=False,
            measure_wilson_loop_mu_nu=False,
            extract_t0=False,
            t0_target=0.1,
            obs_filename="gradient_flow_obs.dat",
            W_temp_filename="gradient_flow_wtemp.dat",
            W_mu_nu_filename="gradient_flow_w_mu_nu.dat",
            t0_filename="gradient_flow_t0.dat",
        )

    payload = {
        "metro": metro_dict,
        "gauge": asdict(gauge),
        "gradient_flow": asdict(gradient_flow),
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.md5(raw.encode("utf-8")).hexdigest()


def _group_key_for_run(path: str) -> Optional[str]:
    yaml_path = Path(path) / "input.yaml"
    if not yaml_path.exists():
        return None
    try:
        metro, gauge = load_params(str(yaml_path))
        gradient_flow = load_gradient_flow_params(str(yaml_path))
    except Exception:
        return None
    return _group_key_from_params(metro, gauge, gradient_flow)


def _build_group_index(root: Path) -> Dict[str, List[str]]:
    index: Dict[str, List[str]] = {}
    if not root.exists():
        return index

    for dirpath, _, filenames in os.walk(root):
        if "input.yaml" not in filenames:
            continue
        key = _group_key_for_run(dirpath)
        if key is None:
            continue
        index.setdefault(key, []).append(os.path.abspath(dirpath))

    for key in index:
        index[key] = sorted(set(index[key]))
    return index


def _get_group_index() -> Dict[str, List[str]]:
    global _GROUP_INDEX_CACHE, _GROUP_INDEX_ROOT
    if _GROUP_INDEX_CACHE is None or _GROUP_INDEX_ROOT != DATA_ROOT:
        _GROUP_INDEX_CACHE = _build_group_index(DATA_ROOT)
        _GROUP_INDEX_ROOT = DATA_ROOT
    return _GROUP_INDEX_CACHE


def _discover_equivalent_runs(path: str) -> Tuple[str, List[str]]:
    abs_path = os.path.abspath(path)
    group_key = _group_key_for_run(abs_path)
    if group_key is None:
        return f"run_{get_run_id(abs_path)}", [abs_path]

    index = _get_group_index()
    grouped_paths = list(index.get(group_key, []))
    if abs_path not in grouped_paths:
        grouped_paths.append(abs_path)
    grouped_paths = sorted(set(grouped_paths))
    return f"group_{group_key}", grouped_paths


def _resolve_worker_count(requested: Optional[int], task_count: int, default: Optional[int] = None) -> int:
    if task_count <= 1:
        return 1

    if requested is None:
        requested = default if default is not None else (os.cpu_count() or 1)

    return max(1, min(int(requested), task_count))


def _load_single_compact_w_temp(
    run_path: str,
    thermalization_steps: Optional[int] = None,
) -> Tuple[str, Optional[do.CompactWilsonData]]:
    run_dir = Path(run_path)
    candidates = [run_dir / "W_temp.out"]
    try:
        flow_params = load_gradient_flow_params(str(run_dir / "input.yaml"))
        candidates.append(run_dir / flow_params.W_temp_filename)
    except Exception:
        pass

    min_step = THERMALIZATION_STEPS if thermalization_steps is None else int(thermalization_steps)
    compact_parts = [
        do.load_compact_wilson_file(str(path), min_step=min_step)
        for path in candidates
        if path.exists()
    ]
    compact_parts = [part for part in compact_parts if part is not None]
    if not compact_parts:
        return run_path, None

    flow_order: list[do.FlowKey] = []
    wilson_by_flow_pair: dict[do.FlowKey, np.ndarray] = {}
    for compact in compact_parts:
        for key in compact.flow_pair_order:
            if key in wilson_by_flow_pair:
                raise ValueError(f"Duplicate Wilson-loop key {key} while loading {run_path}")
            flow_order.append(key)
            wilson_by_flow_pair[key] = compact.wilson_by_flow_pair[key]

    return run_path, do.CompactWilsonData(
        f"{run_path}/combined_w_temp",
        flow_pair_order=flow_order,
        wilson_by_flow_pair=wilson_by_flow_pair,
    )


def _load_single_gradient_flow_obs(
    run_path: str,
    thermalization_steps: Optional[int] = None,
) -> Tuple[str, Optional[Dict[float, Dict[str, np.ndarray]]], Optional[float]]:
    run_dir = Path(run_path)
    try:
        flow_params = load_gradient_flow_params(str(run_dir / "input.yaml"))
    except Exception:
        return run_path, None, None

    obs_path = run_dir / flow_params.obs_filename
    if not obs_path.exists():
        return run_path, None, float(flow_params.t0_target)

    fd = do.FileData(str(obs_path))
    fd.read_file()
    fd.align_lengths()
    min_step = THERMALIZATION_STEPS if thermalization_steps is None else int(thermalization_steps)
    if min_step > 0:
        fd.remove_thermalization(min_step)
    fd.align_lengths()

    if not fd.observables:
        return run_path, None, float(flow_params.t0_target)

    try:
        t_values = np.asarray(fd.get("t_over_a2").values, dtype=float)
        ehat = np.asarray(fd.get("Ehat_clover").values, dtype=float)
        t2e = np.asarray(fd.get("t2E_clover").values, dtype=float)
    except ValueError:
        return run_path, None, float(flow_params.t0_target)

    grouped: Dict[float, Dict[str, list[float]]] = {}
    for t_val, e_val, t2e_val in zip(t_values, ehat, t2e):
        key = float(round(float(t_val), 12))
        row = grouped.setdefault(key, {"Ehat_clover": [], "t2E_clover": []})
        row["Ehat_clover"].append(float(e_val))
        row["t2E_clover"].append(float(t2e_val))

    flow_data = {
        t_val: {
            "Ehat_clover": np.asarray(values["Ehat_clover"], dtype=np.float32),
            "t2E_clover": np.asarray(values["t2E_clover"], dtype=np.float32),
        }
        for t_val, values in grouped.items()
    }
    return run_path, flow_data, float(flow_params.t0_target)


def _load_combined_gradient_flow_obs(
    run_paths: List[str],
    load_workers: int = 1,
    thermalization_steps: Optional[int] = None,
    thermalization_steps_by_run: Optional[Dict[str, int]] = None,
) -> Tuple[Optional[Dict[float, Dict[str, np.ndarray]]], Dict[str, Any]]:
    cut_by_run = {
        os.path.abspath(path): int(value)
        for path, value in (thermalization_steps_by_run or {}).items()
    }

    def cut_for_run(run_path: str) -> Optional[int]:
        return cut_by_run.get(os.path.abspath(run_path), thermalization_steps)

    runs_with_flow = 0
    target = None
    common_times: set[float] | None = None
    chunks: Dict[float, Dict[str, list[np.ndarray]]] = {}

    for run_path in run_paths:
        _, flow_data, run_target = _load_single_gradient_flow_obs(
            run_path,
            thermalization_steps=cut_for_run(run_path),
        )
        if run_target is not None:
            target = run_target if target is None else target
        if not flow_data:
            continue
        runs_with_flow += 1
        current_times = set(flow_data)
        if common_times is None:
            common_times = set(current_times)
        else:
            common_times &= current_times
        for t_val in current_times:
            slot = chunks.setdefault(t_val, {"Ehat_clover": [], "t2E_clover": []})
            slot["Ehat_clover"].append(flow_data[t_val]["Ehat_clover"])
            slot["t2E_clover"].append(flow_data[t_val]["t2E_clover"])

    if not common_times:
        return None, {
            "n_runs_with_gradient_flow": runs_with_flow,
            "t0_target": target,
        }

    combined = {
        t_val: {
            name: np.concatenate(chunks[t_val][name]).astype(np.float32, copy=False)
            for name in ("Ehat_clover", "t2E_clover")
        }
        for t_val in sorted(common_times)
    }
    first = next(iter(combined.values()))
    return combined, {
        "n_runs_with_gradient_flow": runs_with_flow,
        "available_flow_times": sorted(float(t) for t in common_times),
        "n_configurations_after_cut": int(len(first["t2E_clover"])),
        "t0_target": target,
    }


def _bootstrap_series_matrix(
    series_by_row: list[np.ndarray],
    block_size: int,
    n_bootstrap: int,
    seed: int = 42,
) -> np.ndarray:
    if not series_by_row:
        return np.empty((0, n_bootstrap), dtype=float)
    n_samples = min(len(row) for row in series_by_row)
    if n_samples <= 0:
        return np.empty((len(series_by_row), n_bootstrap), dtype=float)
    if any(len(row) != n_samples for row in series_by_row):
        series_by_row = [np.asarray(row, dtype=float)[:n_samples] for row in series_by_row]

    block_size = max(1, int(block_size))
    starts = np.arange(0, n_samples, block_size, dtype=np.int64)
    lengths = np.minimum(block_size, n_samples - starts).astype(np.float64)
    sums = np.empty((len(starts), len(series_by_row)), dtype=np.float64)
    for idx, row in enumerate(series_by_row):
        sums[:, idx] = np.add.reduceat(np.asarray(row, dtype=np.float64), starts)

    rng = np.random.default_rng(seed)
    n_blocks = len(starts)
    out = np.empty((len(series_by_row), n_bootstrap), dtype=float)
    for boot_idx in range(n_bootstrap):
        sampled = rng.integers(0, n_blocks, size=n_blocks)
        out[:, boot_idx] = sums[sampled].sum(axis=0) / lengths[sampled].sum()
    return out


def _trim_series_to_common_length(series_by_row: list[np.ndarray]) -> tuple[list[np.ndarray], int, list[int]]:
    lengths = [int(len(row)) for row in series_by_row]
    if not lengths:
        return [], 0, []
    common_length = min(lengths)
    if common_length <= 0:
        return [np.asarray(row, dtype=float)[:0] for row in series_by_row], 0, lengths
    return [np.asarray(row, dtype=float)[:common_length] for row in series_by_row], common_length, lengths


def _interpolate_t0(flow_times: np.ndarray, t2e_values: np.ndarray, target: float) -> float:
    valid = np.isfinite(flow_times) & np.isfinite(t2e_values)
    x = flow_times[valid]
    y = t2e_values[valid]
    if len(x) < 2:
        return np.nan
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    shifted = y - float(target)
    for idx in range(len(x) - 1):
        y1 = shifted[idx]
        y2 = shifted[idx + 1]
        if y1 == 0:
            return float(x[idx])
        if y1 * y2 <= 0 and y2 != y1:
            return float(x[idx] + (-y1) * (x[idx + 1] - x[idx]) / (y2 - y1))
    if shifted[-1] == 0:
        return float(x[-1])
    return np.nan


def _interpolate_nearest_target_crossing(flow_times: np.ndarray, t2e_values: np.ndarray, target: float) -> float:
    valid = np.isfinite(flow_times) & np.isfinite(t2e_values)
    x = np.asarray(flow_times, dtype=float)[valid]
    y = np.asarray(t2e_values, dtype=float)[valid]
    if len(x) < 1:
        return np.nan

    exact = np.flatnonzero(y == float(target))
    if exact.size > 0:
        return float(x[exact[0]])

    below = np.flatnonzero(y < float(target))
    above = np.flatnonzero(y > float(target))
    if below.size == 0 or above.size == 0:
        return np.nan

    below_idx = below[np.argmin(np.abs(y[below] - float(target)))]
    above_idx = above[np.argmin(np.abs(y[above] - float(target)))]
    y1 = y[below_idx]
    y2 = y[above_idx]
    if y2 == y1:
        return np.nan
    return float(x[below_idx] + (float(target) - y1) * (x[above_idx] - x[below_idx]) / (y2 - y1))


def _std_finite(values: np.ndarray) -> float:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    return float(np.std(finite)) if finite.size > 0 else np.nan


def _weighted_linear_crossing_for_indices(
    flow_times: np.ndarray,
    t2e_values: np.ndarray,
    t2e_errors: np.ndarray,
    target: float,
    indices: np.ndarray,
) -> float:
    idx = np.asarray(indices, dtype=int)
    if idx.size < 2:
        return np.nan

    x = np.asarray(flow_times, dtype=float)[idx]
    y = np.asarray(t2e_values, dtype=float)[idx]
    err = np.asarray(t2e_errors, dtype=float)[idx]
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]
    err = err[valid]
    if x.size < 2:
        return np.nan

    weights = None
    if err.size == x.size and np.all(np.isfinite(err)) and np.all(err > 0):
        weights = 1.0 / err
    try:
        slope, intercept = np.polyfit(x, y, deg=1, w=weights)
    except (ValueError, np.linalg.LinAlgError):
        return np.nan
    if not np.isfinite(slope) or abs(float(slope)) <= np.finfo(float).eps:
        return np.nan
    return float((float(target) - float(intercept)) / float(slope))


def _weighted_linear_crossing(
    flow_times: np.ndarray,
    t2e_values: np.ndarray,
    t2e_errors: np.ndarray,
    target: float,
    points_per_side: int = 2,
) -> tuple[float, np.ndarray]:
    y = np.asarray(t2e_values, dtype=float)
    valid = np.isfinite(flow_times) & np.isfinite(y)
    below = np.flatnonzero(valid & (y < float(target)))
    above = np.flatnonzero(valid & (y > float(target)))
    if below.size == 0 or above.size == 0:
        return np.nan, np.asarray([], dtype=int)

    below = below[np.argsort(np.abs(y[below] - float(target)))[:points_per_side]]
    above = above[np.argsort(np.abs(y[above] - float(target)))[:points_per_side]]
    indices = np.asarray(sorted(np.concatenate([below, above]).tolist()), dtype=int)
    return (
        _weighted_linear_crossing_for_indices(flow_times, y, t2e_errors, target, indices),
        indices,
    )


def summarize_gradient_flow_obs(
    flow_data: Optional[Dict[float, Dict[str, np.ndarray]]],
    *,
    t0_target: float | None,
    block_size: int,
    n_bootstrap: int,
    include_bootstrap: bool = False,
) -> Dict[str, Any]:
    if not flow_data:
        return {}
    times = np.asarray(sorted(flow_data), dtype=float)
    e_series = [np.asarray(flow_data[float(t)]["Ehat_clover"], dtype=float) for t in times]
    t2e_series = [np.asarray(flow_data[float(t)]["t2E_clover"], dtype=float) for t in times]
    e_series, n_common_e, e_lengths = _trim_series_to_common_length(e_series)
    t2e_series, n_common_t2e, t2e_lengths = _trim_series_to_common_length(t2e_series)
    n_common = min(n_common_e, n_common_t2e)
    if n_common <= 0:
        return {}
    if n_common_e != n_common:
        e_series = [row[:n_common] for row in e_series]
    if n_common_t2e != n_common:
        t2e_series = [row[:n_common] for row in t2e_series]
    e_means = np.asarray([np.mean(row) for row in e_series], dtype=float)
    t2e_means = np.asarray([np.mean(row) for row in t2e_series], dtype=float)
    e_boot = _bootstrap_series_matrix(e_series, block_size, n_bootstrap)
    t2e_boot = _bootstrap_series_matrix(t2e_series, block_size, n_bootstrap)
    e_errs = np.asarray([_std_finite(e_boot[idx]) for idx in range(e_boot.shape[0])], dtype=float)
    t2e_errs = np.asarray([_std_finite(t2e_boot[idx]) for idx in range(t2e_boot.shape[0])], dtype=float)
    target = 0.1 if t0_target is None else float(t0_target)
    t0 = _interpolate_t0(times, t2e_means, target)
    t0_boot = np.asarray([
        _interpolate_t0(times, t2e_boot[:, idx], target)
        for idx in range(t2e_boot.shape[1])
    ], dtype=float)
    finite_t0 = t0_boot[np.isfinite(t0_boot)]
    fixed_target = 0.1
    t2e_0p1 = _interpolate_nearest_target_crossing(times, t2e_means, fixed_target)
    t2e_0p1_boot = np.asarray([
        _interpolate_nearest_target_crossing(times, t2e_boot[:, idx], fixed_target)
        for idx in range(t2e_boot.shape[1])
    ], dtype=float)
    finite_t2e_0p1 = t2e_0p1_boot[np.isfinite(t2e_0p1_boot)]
    t2e_0p1_fit, fit_indices = _weighted_linear_crossing(
        times,
        t2e_means,
        t2e_errs,
        fixed_target,
    )
    t2e_0p1_fit_boot = np.asarray([
        _weighted_linear_crossing_for_indices(times, t2e_boot[:, idx], t2e_errs, fixed_target, fit_indices)
        for idx in range(t2e_boot.shape[1])
    ], dtype=float)
    finite_t2e_0p1_fit = t2e_0p1_fit_boot[np.isfinite(t2e_0p1_fit_boot)]

    payload: Dict[str, Any] = {
        "available_flow_times": [float(t) for t in times],
        "n_configurations_used": int(n_common),
        "n_configurations_by_flow_time": {
            f"{float(t):.12g}": {
                "Ehat_clover": int(e_len),
                "t2E_clover": int(t2e_len),
            }
            for t, e_len, t2e_len in zip(times, e_lengths, t2e_lengths)
        },
        "truncated_to_common_length": bool(
            any(length != n_common for length in e_lengths)
            or any(length != n_common for length in t2e_lengths)
        ),
        "Ehat_clover": {f"{float(t):.12g}": float(v) for t, v in zip(times, e_means)},
        "Ehat_clover_err": {
            f"{float(t):.12g}": float(e_errs[idx]) if np.isfinite(e_errs[idx]) else None
            for idx, t in enumerate(times)
        },
        "t2E_clover": {f"{float(t):.12g}": float(v) for t, v in zip(times, t2e_means)},
        "t2E_clover_err": {
            f"{float(t):.12g}": float(t2e_errs[idx]) if np.isfinite(t2e_errs[idx]) else None
            for idx, t in enumerate(times)
        },
        "t0_target": target,
        "t0": float(t0) if np.isfinite(t0) else None,
        "t0_err": float(np.std(finite_t0)) if finite_t0.size > 0 else None,
        "t_over_a2_at_t2E_clover_0p1": float(t2e_0p1) if np.isfinite(t2e_0p1) else None,
        "t_over_a2_at_t2E_clover_0p1_err": (
            float(np.std(finite_t2e_0p1)) if finite_t2e_0p1.size > 0 else None
        ),
        "t_over_a2_at_t2E_clover_0p1_weighted_fit": (
            float(t2e_0p1_fit) if np.isfinite(t2e_0p1_fit) else None
        ),
        "t_over_a2_at_t2E_clover_0p1_weighted_fit_err": (
            float(np.std(finite_t2e_0p1_fit)) if finite_t2e_0p1_fit.size > 0 else None
        ),
        "t_over_a2_at_t2E_clover_0p1_weighted_fit_points": [
            {
                "t_over_a2": float(times[idx]),
                "t2E_clover": float(t2e_means[idx]),
                "t2E_clover_err": float(t2e_errs[idx]) if np.isfinite(t2e_errs[idx]) else None,
            }
            for idx in fit_indices
        ],
    }
    if include_bootstrap:
        payload["bootstrap_samples"] = {
            "Ehat_clover": e_boot,
            "t2E_clover": t2e_boot,
            "t0": t0_boot,
            "t_over_a2_at_t2E_clover_0p1": t2e_0p1_boot,
            "t_over_a2_at_t2E_clover_0p1_weighted_fit": t2e_0p1_fit_boot,
        }
    return payload


def _iter_loaded_compact_w_temp(
    run_paths: List[str],
    load_workers: int,
    thermalization_steps: Optional[int] = None,
    thermalization_steps_by_run: Optional[Dict[str, int]] = None,
):
    cut_by_run = {
        os.path.abspath(path): int(value)
        for path, value in (thermalization_steps_by_run or {}).items()
    }

    def cut_for_run(run_path: str) -> Optional[int]:
        return cut_by_run.get(os.path.abspath(run_path), thermalization_steps)

    if load_workers <= 1:
        for run_path in run_paths:
            yield _load_single_compact_w_temp(run_path, thermalization_steps=cut_for_run(run_path))
        return

    with concurrent.futures.ThreadPoolExecutor(max_workers=load_workers) as executor:
        future_to_index: Dict[concurrent.futures.Future, int] = {}
        ready_results: Dict[int, Tuple[str, Optional[do.CompactWilsonData]]] = {}
        next_submit = 0
        next_yield = 0

        def submit_one(index: int):
            future = executor.submit(
                _load_single_compact_w_temp,
                run_paths[index],
                cut_for_run(run_paths[index]),
            )
            future_to_index[future] = index

        initial = min(load_workers, len(run_paths))
        for _ in range(initial):
            submit_one(next_submit)
            next_submit += 1

        while next_yield < len(run_paths):
            if next_yield in ready_results:
                yield ready_results.pop(next_yield)
                if next_submit < len(run_paths):
                    submit_one(next_submit)
                    next_submit += 1
                next_yield += 1
                continue

            done, _ = concurrent.futures.wait(
                tuple(future_to_index),
                return_when=concurrent.futures.FIRST_COMPLETED,
            )
            for future in done:
                index = future_to_index.pop(future)
                ready_results[index] = future.result()


def _calculate_w_rt_variables(
    file_data: do.FileData,
    block_size: int,
    available_pairs: List[Tuple[int, int]],
    calc_workers: Optional[int],
    verbose: bool = False,
    prefix: str = "",
) -> Tuple[Dict[str, float], Dict[str, float], Dict[Any, do.VariableData]]:
    def vprint(msg: str):
        if verbose:
            import sys
            print(f"{prefix} {msg}", file=sys.stderr, flush=True)

    pair_list = [(int(r), int(t)) for r, t in available_pairs]
    worker_count = _resolve_worker_count(calc_workers, len(pair_list))
    w_values: Dict[str, float] = {}
    w_errors: Dict[str, float] = {}
    w_cache: Dict[Any, do.VariableData] = {}

    if not pair_list:
        return w_values, w_errors, w_cache

    vprint(f"Calculating {len(pair_list)} W(R,T) value(s) with {worker_count} worker(s)...")

    thread_state = threading.local()

    def get_calc() -> Calculator:
        calc = getattr(thread_state, "calc", None)
        if calc is None:
            calc = Calculator(file_data, n_bootstrap=DEFAULT_N_BOOTSTRAP, step_size=block_size)
            thread_state.calc = calc
        return calc

    def task(pair: Tuple[int, int]) -> Tuple[Tuple[int, int], Optional[do.VariableData]]:
        r_val, t_val = pair
        try:
            var = get_calc().get_variable("W_R_T", R=r_val, T=t_val)
            return pair, var
        except Exception:
            return pair, None

    if worker_count == 1:
        iterator = map(task, pair_list)
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=worker_count) as executor:
            iterator = executor.map(task, pair_list)
            for pair, var in iterator:
                if var is None:
                    continue
                r_val, t_val = pair
                w_cache[make_key("W_R_T", {"R": r_val, "T": t_val})] = var
                w_val = var.get()
                if w_val is not None:
                    w_values[f"{r_val},{t_val}"] = float(w_val)
                    if var.err() is not None:
                        w_errors[f"{r_val},{t_val}"] = float(var.err())
            return w_values, w_errors, w_cache

    for pair, var in iterator:
        if var is None:
            continue
        r_val, t_val = pair
        w_cache[make_key("W_R_T", {"R": r_val, "T": t_val})] = var
        w_val = var.get()
        if w_val is not None:
            w_values[f"{r_val},{t_val}"] = float(w_val)
            if var.err() is not None:
                w_errors[f"{r_val},{t_val}"] = float(var.err())

    return w_values, w_errors, w_cache


def _calculate_v_r_variables(
    file_data: do.FileData,
    block_size: int,
    unique_Rs: List[int],
    calc_workers: Optional[int],
    seed_cache: Optional[Dict[Any, do.VariableData]] = None,
    extra_params: Optional[Dict[str, Any]] = None,
    verbose: bool = False,
    prefix: str = "",
) -> Tuple[Dict[str, float], Dict[str, float], Dict[Any, do.VariableData]]:
    def vprint(msg: str):
        if verbose:
            import sys
            print(f"{prefix} {msg}", file=sys.stderr, flush=True)

    r_values = [int(r) for r in unique_Rs]
    worker_count = _resolve_worker_count(calc_workers, len(r_values))
    v_params = dict(extra_params or {})
    shared_cache = dict(seed_cache or {})
    potentials: Dict[str, float] = {}
    potential_errors: Dict[str, float] = {}
    v_cache: Dict[Any, do.VariableData] = {}

    if not r_values:
        return potentials, potential_errors, v_cache

    vprint(f"Calculating {len(r_values)} V(R) value(s) with {worker_count} worker(s)...")

    thread_state = threading.local()

    def build_calc() -> Calculator:
        calc = Calculator(file_data, n_bootstrap=DEFAULT_N_BOOTSTRAP, step_size=block_size)
        if shared_cache:
            calc.variables.update(shared_cache)
        return calc

    def get_calc() -> Calculator:
        calc = getattr(thread_state, "calc", None)
        if calc is None:
            calc = build_calc()
            thread_state.calc = calc
        return calc

    def task(r_val: int) -> Tuple[int, Optional[do.VariableData]]:
        try:
            var = get_calc().get_variable("V_R", R=r_val, **v_params)
            return r_val, var
        except Exception:
            return r_val, None

    if worker_count == 1:
        iterator = map(task, r_values)
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=worker_count) as executor:
            iterator = executor.map(task, r_values)
            for r_val, var in iterator:
                if var is None:
                    continue
                v_cache[make_key("V_R", {"R": r_val, **v_params})] = var
                val = var.get()
                if val is not None and not np.isnan(val):
                    potentials[str(r_val)] = float(val)
                    potential_errors[str(r_val)] = float(var.err())
            return potentials, potential_errors, v_cache

    for r_val, var in iterator:
        if var is None:
            continue
        v_cache[make_key("V_R", {"R": r_val, **v_params})] = var
        val = var.get()
        if val is not None and not np.isnan(val):
            potentials[str(r_val)] = float(val)
            potential_errors[str(r_val)] = float(var.err())

    return potentials, potential_errors, v_cache


def _load_combined_w_temp(
    run_paths: List[str],
    verbose: bool = False,
    prefix: str = "",
    load_workers: int = 1,
    thermalization_steps: Optional[int] = None,
    thermalization_steps_by_run: Optional[Dict[str, int]] = None,
) -> Tuple[Optional[do.FileData], Dict[str, Any]]:
    def vprint(msg: str):
        if verbose:
            import sys
            print(f"{prefix} {msg}", file=sys.stderr, flush=True)

    worker_count = _resolve_worker_count(load_workers, len(run_paths), default=1)
    runs_with_w_temp = 0
    resolved_cuts_by_run = {
        os.path.abspath(path): (
            int(thermalization_steps_by_run[os.path.abspath(path)])
            if thermalization_steps_by_run is not None and os.path.abspath(path) in thermalization_steps_by_run
            else (THERMALIZATION_STEPS if thermalization_steps is None else int(thermalization_steps))
        )
        for path in run_paths
    }
    unique_cuts = sorted(set(resolved_cuts_by_run.values()))

    def iter_compact_files():
        nonlocal runs_with_w_temp
        for run_path, compact in _iter_loaded_compact_w_temp(
            run_paths,
            worker_count,
            thermalization_steps=thermalization_steps,
            thermalization_steps_by_run=resolved_cuts_by_run,
        ):
            if compact is None:
                continue
            runs_with_w_temp += 1
            vprint(f"Loaded data from {Path(run_path).name}...")
            yield compact

    if worker_count > 1:
        vprint(f"Loading up to {worker_count} W_temp file(s) in parallel...")
    vprint("Combining W_temp data incrementally...")
    vprint(f"Total runs in group: {len(run_paths)}")
    combined = do.combine_compact_wilson_data(
        iter_compact_files(),
        source_name="W_temp_combined",
    )
    vprint("Combination complete.")

    metadata = {
        "n_runs_in_group": len(run_paths),
        "n_runs_with_w_temp": runs_with_w_temp,
        "n_w_temp_files": runs_with_w_temp,
        "n_samples_after_cut": 0,
        "n_configurations_after_cut": 0,
        "thermalization_steps": unique_cuts[0] if len(unique_cuts) == 1 else None,
        "thermalization_steps_by_run": resolved_cuts_by_run,
        "load_workers": worker_count,
    }

    if combined is not None:
        pair_sample_counts = getattr(combined, "pair_sample_counts", {})
        sample_lengths = [int(value) for value in pair_sample_counts.values()]
        n_configurations_after_cut = int(max(sample_lengths, default=getattr(combined, "n_configurations", 0)))

        def format_flow_key(flow_time: float | None, r_val: int, t_val: int) -> str:
            flow_text = "none" if flow_time is None else f"{float(flow_time):g}"
            return f"{flow_text},{int(r_val)},{int(t_val)}"

        metadata["n_configurations_after_cut"] = n_configurations_after_cut
        metadata["min_configurations_per_wilson_loop"] = int(min(sample_lengths, default=0))
        metadata["max_configurations_per_wilson_loop"] = int(max(sample_lengths, default=0))
        metadata["n_samples_after_cut"] = int(sum(sample_lengths))
        metadata["wilson_loop_sample_counts"] = {
            format_flow_key(flow_time, r_val, t_val): int(count)
            for (flow_time, r_val, t_val), count in sorted(
                pair_sample_counts.items(),
                key=lambda item: (-1.0 if item[0][0] is None else float(item[0][0]), item[0][1], item[0][2]),
            )
        }
        metadata["available_wilson_flow_times"] = [
            None if value is None else float(value)
            for value in combined.available_flow_times()
        ]

    return combined, metadata


def evaluate_run(
    file_data: do.FileData,
    input_dir: Path,
    sommer_target: float = DEFAULT_SOMMER_TARGET,
    verbose: bool = False,
    prefix: str = "",
    calc_workers: Optional[int] = None,
    analysis_options: Optional[Dict[str, Any]] = None,
    gradient_flow_data: Optional[Dict[float, Dict[str, np.ndarray]]] = None,
    gradient_flow_metadata: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    def vprint(msg: str):
        if verbose:
            import sys
            print(f"{prefix} {msg}", file=sys.stderr, flush=True)

    analysis_options = _resolve_analysis_options(analysis_options)
    r0_t_min = analysis_options["r0_t_min"]

    if (
        not file_data.observables
        and not getattr(file_data, "wilson_by_pair", None)
        and not getattr(file_data, "wilson_by_flow_pair", None)
    ):
        return {"error": "No combined W_temp observables available"}

    vprint("Calculating tau_int...")
    tau = 0.5
    if any(o.name in {"plaquette", "retrace"} for o in file_data.observables):
        calc_pre = Calculator(file_data, n_bootstrap=DEFAULT_N_BOOTSTRAP)
        try:
            tau_var = calc_pre.get_variable("tau_int", obs_name="plaquette")
            tau = tau_var.get()
        except KeyError:
            tau = 0.5

    block_size = max(1, int(np.ceil(2 * tau)))
    block_size = max(block_size, DEFAULT_BLOCK_SIZE)
    vprint(f"Setting up Calculator with block_size={block_size}...")
    calc = Calculator(file_data, n_bootstrap=DEFAULT_N_BOOTSTRAP, step_size=block_size)

    vprint("Extracting unique R and T...")
    unique_Rs = calc.get_unique_Rs()
    available_pairs = calc.get_available_pairs()
    unique_Ts = calc.get_unique_Ts()
    analysis_settings = _build_analysis_settings(
        unique_Rs=unique_Rs,
        unique_Ts=unique_Ts,
        tau=tau,
        block_size=block_size,
        sommer_target=sommer_target,
        r0_t_min=r0_t_min,
    )
    flow_summary = summarize_gradient_flow_obs(
        gradient_flow_data,
        t0_target=(gradient_flow_metadata or {}).get("t0_target"),
        block_size=block_size,
        n_bootstrap=DEFAULT_N_BOOTSTRAP,
        include_bootstrap=False,
    )
    if flow_summary:
        analysis_settings["gradient_flow"] = {
            "available_flow_times": flow_summary.get("available_flow_times", []),
            "t0_target": flow_summary.get("t0_target"),
        }

    all_w, all_w_err, w_cache = _calculate_w_rt_variables(
        file_data,
        block_size,
        available_pairs,
        calc_workers,
        verbose=verbose,
        prefix=prefix,
    )
    calc.variables.update(w_cache)

    potentials, potential_errors, potential_cache = _calculate_v_r_variables(
        file_data,
        block_size,
        unique_Rs,
        calc_workers,
        seed_cache=dict(calc.variables),
        verbose=verbose,
        prefix=prefix,
    )
    calc.variables.update(potential_cache)

    vprint("Calculating r0...")
    r0 = None
    r0_err = None
    cornell_params = None
    try:
        r0_fit_Rs = [int(r) for r in unique_Rs if int(r) >= 2]
        _, _, r0_fit_cache = _calculate_v_r_variables(
            file_data,
            block_size,
            r0_fit_Rs,
            calc_workers,
            seed_cache=dict(calc.variables),
            extra_params={"t_min": r0_t_min, "t_max": DEFAULT_R0_T_MAX},
            verbose=verbose,
            prefix=prefix,
        )
        calc.variables.update(r0_fit_cache)
        r0_var = calc.get_variable(
            "r0",
            t_min=r0_t_min,
            t_max=DEFAULT_R0_T_MAX,
            target_force=sommer_target,
            r_min=DEFAULT_R0_R_MIN,
        )
        r0 = r0_var.get()
        if np.isnan(r0):
            r0 = None
        else:
            r0_err = r0_var.err()
            cornell_params = r0_var.parameters.get("cornell_params", None)
    except Exception:
        pass

    lattice_spacing = None
    a_err = None
    if r0 is not None:
        lattice_spacing = 0.5 / r0
        a_err = (0.5 / (r0**2)) * r0_err if r0_err else 0.0

    epsilon_bar = None
    epsilon_bar_err = None
    metro = None
    try:
        yaml_path = input_dir / "input.yaml"
        metro, _ = load_params(str(yaml_path))
        if r0 is not None:
            epsilon_bar_var = calc.get_variable(
                "epsilon_bar",
                epsilon1=metro.epsilon1,
                beta=metro.beta,
                t_min=r0_t_min,
                t_max=DEFAULT_R0_T_MAX,
                target_force=sommer_target,
                r_min=DEFAULT_R0_R_MIN,
            )
            if epsilon_bar_var.get() is not None and not np.isnan(epsilon_bar_var.get()):
                epsilon_bar = epsilon_bar_var.get()
                epsilon_bar_err = epsilon_bar_var.err()
    except Exception:
        pass

    vprint("Calculating chi and F_chi...")
    all_chi = {}
    all_chi_err = {}
    all_f_chi = {}
    all_f_chi_err = {}
    r0_chi = None
    r0_chi_err = None
    creutz_status = "ok"
    chi_t_large = DEFAULT_R0_CHI_T_LARGE

    try:
        chi_pairs = []
        for r in unique_Rs:
            if r + 1 not in unique_Rs:
                continue
            for t in unique_Ts:
                if t + 1 not in unique_Ts:
                    continue
                chi_pairs.append((r + 0.5, t + 0.5))

        if not chi_pairs:
            creutz_status = "not enough adjacent L values for standard Creutz ratios"

        for r_p, t_p in chi_pairs:
            try:
                chi_var = calc.get_variable("chi", R=r_p, T=t_p)
                chi_val = chi_var.get()
                if chi_val is not None and not np.isnan(chi_val):
                    key = f"{float(r_p):g},{float(t_p):g}"
                    all_chi[key] = float(chi_val)
                    chi_err = chi_var.err()
                    if chi_err is not None and np.isfinite(chi_err):
                        all_chi_err[key] = float(chi_err)
            except Exception:
                continue

        r_for_force = sorted({p[0] for p in chi_pairs if p[1] == chi_t_large + 0.5})
        for r_p in r_for_force:
            try:
                f_var = calc.get_variable("F_chi", R=r_p, t_large=chi_t_large)
                f_val = f_var.get()
                if f_val is not None and not np.isnan(f_val):
                    key = f"{float(r_p):g}"
                    all_f_chi[key] = float(f_val)
                    f_err = f_var.err()
                    if f_err is not None and np.isfinite(f_err):
                        all_f_chi_err[key] = float(f_err)
            except Exception:
                continue

        r0_chi_var = calc.get_variable(
            "r0_chi",
            t_large=chi_t_large,
            target_force=sommer_target,
            max_rel_err=DEFAULT_R0_CHI_MAX_REL_ERR,
            use_weighted_fit=DEFAULT_R0_CHI_USE_WEIGHTED_FIT,
            fit_window=DEFAULT_R0_CHI_FIT_WINDOW,
            discard_negative=DEFAULT_R0_CHI_DISCARD_NEGATIVE,
            r_min=DEFAULT_R0_CHI_R_MIN,
        )
        r0_chi = r0_chi_var.get()
        if np.isnan(r0_chi):
            r0_chi = None
        else:
            r0_chi_err = r0_chi_var.err()
    except Exception:
        if creutz_status == "ok":
            creutz_status = "failed"

    vprint("Calculating volume_r0 and length...")
    volume_r0 = None
    volume_r0_err = None
    length = None
    length_err = None
    try:
        if metro is not None and r0 is not None:
            vol_var = calc.get_variable(
                "volume_r0",
                L0=metro.L0,
                L1=metro.L1,
                L2=metro.L2,
                L3=metro.L3,
                t_min=r0_t_min,
                t_max=DEFAULT_R0_T_MAX,
                target_force=sommer_target,
                r_min=DEFAULT_R0_R_MIN,
            )
            if vol_var.get() is not None and not np.isnan(vol_var.get()):
                volume_r0 = vol_var.get()
                volume_r0_err = vol_var.err()

            len_var = calc.get_variable(
                "length",
                L0=metro.L0,
                t_min=r0_t_min,
                t_max=DEFAULT_R0_T_MAX,
                target_force=sommer_target,
                r_min=DEFAULT_R0_R_MIN,
            )
            if len_var.get() is not None and not np.isnan(len_var.get()):
                length = len_var.get()
                length_err = len_var.err()
    except Exception:
        pass

    vprint("Calculating creutz_P and a_creutz...")
    all_creutz_P = {}
    all_creutz_P_err = {}
    all_a_creutz = {}
    all_a_creutz_err = {}

    try:
        even_Rs = [r for r in unique_Rs if r % 2 == 0 and r > 0]
        for r in even_Rs:
            try:
                p_var = calc.get_variable("creutz_P", R=int(r))
                p_val = p_var.get()
                if p_val is not None and not np.isnan(p_val):
                    all_creutz_P[str(int(r))] = float(p_val)
                    if p_var.err() is not None:
                        all_creutz_P_err[str(int(r))] = float(p_var.err())

                a_var = calc.get_variable("a_creutz", R=int(r))
                a_val = a_var.get()
                if a_val is not None and not np.isnan(a_val):
                    all_a_creutz[str(int(r))] = float(a_val)
                    if a_var.err() is not None:
                        all_a_creutz_err[str(int(r))] = float(a_var.err())
            except Exception:
                continue
    except Exception:
        pass

    return {
        "r0": float(r0) if r0 is not None else None,
        "r0_err": float(r0_err) if r0_err is not None else None,
        "volume_r0": float(volume_r0) if volume_r0 is not None else None,
        "volume_r0_err": float(volume_r0_err) if volume_r0_err is not None else None,
        "length": float(length) if length is not None else None,
        "length_err": float(length_err) if length_err is not None else None,
        "a": float(lattice_spacing) if lattice_spacing is not None else None,
        "a_err": float(a_err) if a_err is not None else None,
        "epsilon_bar": float(epsilon_bar) if epsilon_bar is not None else None,
        "epsilon_bar_err": float(epsilon_bar_err) if epsilon_bar_err is not None else None,
        "tau_int": float(tau),
        "block_size": int(block_size),
        "analysis_settings": analysis_settings,
        "V_R": potentials,
        "V_R_err": potential_errors,
        "W_R_T": all_w,
        "W_R_T_err": all_w_err,
        "r0_chi": float(r0_chi) if r0_chi is not None else None,
        "r0_chi_err": float(r0_chi_err) if r0_chi_err is not None else None,
        "chi": all_chi,
        "chi_err": all_chi_err,
        "F_chi": all_f_chi,
        "F_chi_err": all_f_chi_err,
        "creutz_P": all_creutz_P,
        "creutz_P_err": all_creutz_P_err,
        "a_creutz": all_a_creutz,
        "a_creutz_err": all_a_creutz_err,
        "creutz_status": creutz_status,
        "gradient_flow": flow_summary,
        "gradient_flow_metadata": gradient_flow_metadata or {},
        "plot_meta": {
            "potentials": potentials,
            "potential_errors": potential_errors,
            "cornell_params": cornell_params,
            "chi": all_chi,
            "chi_err": all_chi_err,
            "F_chi": all_f_chi,
            "F_chi_err": all_f_chi_err,
            "gradient_flow": flow_summary,
        },
    }


def get_or_calculate(
    path: str,
    force_recalc: bool = False,
    verbose: bool = True,
    calc_workers: Optional[int] = None,
    load_workers: int = 1,
    combine_equivalent_runs: bool = True,
    analysis_options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    path = os.path.abspath(path)
    run_id = get_run_id(path)
    analysis_options = _resolve_analysis_options(analysis_options)

    prefix = f"[{os.path.basename(path)}]"
    def vprint(msg: str):
        if verbose:
            import sys
            print(f"{prefix} {msg}", file=sys.stderr, flush=True)

    if not os.path.isdir(path):
        return {"error": "Directory not found"}

    if combine_equivalent_runs:
        base_analysis_id, grouped_paths = _discover_equivalent_runs(path)
    else:
        base_analysis_id = f"run_{run_id}"
        grouped_paths = [path]
    analysis_id = _cache_id_for_analysis(base_analysis_id, analysis_options)

    if not force_recalc:
        cached = load_cached_result(analysis_id)
        if cached:
            vprint("Loading data from cache...")
            result = dict(cached)
            result["path"] = path
            result["run_id"] = run_id
            result["analysis_id"] = analysis_id
            return result

    try:
        combined_w_temp, aggregation = _load_combined_w_temp(
            grouped_paths,
            verbose=verbose,
            prefix=prefix,
            load_workers=load_workers,
        )
        if combined_w_temp is None:
            return {"error": "No W_temp data found in equivalent run group"}

        gradient_flow_data, gradient_flow_metadata = _load_combined_gradient_flow_obs(
            grouped_paths,
            load_workers=load_workers,
        )

        result = evaluate_run(
            combined_w_temp,
            Path(path),
            verbose=verbose,
            prefix=prefix,
            calc_workers=calc_workers,
            analysis_options=analysis_options,
            gradient_flow_data=gradient_flow_data,
            gradient_flow_metadata=gradient_flow_metadata,
        )
        result["path"] = path
        result["run_id"] = run_id
        result["analysis_id"] = analysis_id
        result["aggregation"] = aggregation

        if "error" not in result:
            try:
                save_result(analysis_id, result)
            except OSError as exc:
                vprint(f"Warning: Could not save cache: {exc}")
        return result
    except Exception as e:
        return {"error": str(e)}
