import json
import hashlib
import os
import threading
import concurrent.futures
from dataclasses import asdict
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple

import numpy as np

from load_input_yaml import load_params, MetropolisParams, GaugeObservableParams

import data_organizer as do
from calculator import Calculator, make_key

# --- CONFIGURATION ---
DATA_ROOT = Path("../data").resolve()
CALC_RESULT_BASE = DATA_ROOT / "calcResult"
THERMALIZATION_STEPS = 500
CALC_VERSION = "4.4"
GROUP_IGNORE_METRO_FIELDS = {"seed", "nSweep"}

DEFAULT_N_BOOTSTRAP = 200
DEFAULT_SOMMER_TARGET = 1.65
DEFAULT_V_R_T_MIN = 1
DEFAULT_V_R_T_MAX = None
DEFAULT_R0_T_MIN = 5
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
            "t_min": DEFAULT_R0_T_MIN,
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


def _group_key_from_params(metro: MetropolisParams, gauge: GaugeObservableParams) -> str:
    metro_dict = asdict(metro)
    for key in GROUP_IGNORE_METRO_FIELDS:
        metro_dict.pop(key, None)

    payload = {
        "metro": metro_dict,
        "gauge": asdict(gauge),
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.md5(raw.encode("utf-8")).hexdigest()


def _group_key_for_run(path: str) -> Optional[str]:
    yaml_path = Path(path) / "input.yaml"
    if not yaml_path.exists():
        return None
    try:
        metro, gauge = load_params(str(yaml_path))
    except Exception:
        return None
    return _group_key_from_params(metro, gauge)


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


def _load_single_compact_w_temp(run_path: str) -> Tuple[str, Optional[do.CompactWilsonData]]:
    w_temp_path = Path(run_path) / "W_temp.out"
    if not w_temp_path.exists():
        return run_path, None

    compact = do.load_compact_wilson_file(str(w_temp_path), min_step=THERMALIZATION_STEPS)
    return run_path, compact


def _iter_loaded_compact_w_temp(
    run_paths: List[str],
    load_workers: int,
):
    if load_workers <= 1:
        for run_path in run_paths:
            yield _load_single_compact_w_temp(run_path)
        return

    with concurrent.futures.ThreadPoolExecutor(max_workers=load_workers) as executor:
        future_to_index: Dict[concurrent.futures.Future, int] = {}
        ready_results: Dict[int, Tuple[str, Optional[do.CompactWilsonData]]] = {}
        next_submit = 0
        next_yield = 0

        def submit_one(index: int):
            future = executor.submit(_load_single_compact_w_temp, run_paths[index])
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
) -> Tuple[Dict[str, float], Dict[Any, do.VariableData]]:
    def vprint(msg: str):
        if verbose:
            import sys
            print(f"{prefix} {msg}", file=sys.stderr, flush=True)

    pair_list = [(int(r), int(t)) for r, t in available_pairs]
    worker_count = _resolve_worker_count(calc_workers, len(pair_list))
    w_values: Dict[str, float] = {}
    w_cache: Dict[Any, do.VariableData] = {}

    if not pair_list:
        return w_values, w_cache

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
            return w_values, w_cache

    for pair, var in iterator:
        if var is None:
            continue
        r_val, t_val = pair
        w_cache[make_key("W_R_T", {"R": r_val, "T": t_val})] = var
        w_val = var.get()
        if w_val is not None:
            w_values[f"{r_val},{t_val}"] = float(w_val)

    return w_values, w_cache


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
) -> Tuple[Optional[do.FileData], Dict[str, Any]]:
    def vprint(msg: str):
        if verbose:
            import sys
            print(f"{prefix} {msg}", file=sys.stderr, flush=True)

    worker_count = _resolve_worker_count(load_workers, len(run_paths), default=1)
    runs_with_w_temp = 0

    def iter_compact_files():
        nonlocal runs_with_w_temp
        for run_path, compact in _iter_loaded_compact_w_temp(run_paths, worker_count):
            if compact is None:
                continue
            runs_with_w_temp += 1
            vprint(f"Loaded data from {Path(run_path).name}...")
            yield compact

    if worker_count > 1:
        vprint(f"Loading up to {worker_count} W_temp file(s) in parallel...")
    vprint("Combining W_temp data incrementally...")
    combined = do.combine_compact_wilson_data(
        iter_compact_files(),
        source_name="W_temp_combined",
    )

    metadata = {
        "n_runs_in_group": len(run_paths),
        "n_runs_with_w_temp": runs_with_w_temp,
        "n_w_temp_files": runs_with_w_temp,
        "n_samples_after_cut": 0,
        "n_configurations_after_cut": 0,
        "thermalization_steps": THERMALIZATION_STEPS,
        "load_workers": worker_count,
    }

    if combined is not None:
        n_configurations_after_cut = getattr(combined, "n_configurations", 0)
        metadata["n_configurations_after_cut"] = n_configurations_after_cut
        metadata["n_samples_after_cut"] = n_configurations_after_cut * len(getattr(combined, "pair_order", []))

    return combined, metadata


def evaluate_run(
    file_data: do.FileData,
    input_dir: Path,
    sommer_target: float = DEFAULT_SOMMER_TARGET,
    verbose: bool = False,
    prefix: str = "",
    calc_workers: Optional[int] = None,
) -> Dict[str, Any]:
    def vprint(msg: str):
        if verbose:
            import sys
            print(f"{prefix} {msg}", file=sys.stderr, flush=True)

    if not file_data.observables and not getattr(file_data, "wilson_by_pair", None):
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

    vprint(f"Setting up Calculator with block_size={max(1, int(np.ceil(2 * tau)))}...")
    block_size = max(1, int(np.ceil(2 * tau)))
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
    )

    all_w, w_cache = _calculate_w_rt_variables(
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
            extra_params={"t_min": DEFAULT_R0_T_MIN, "t_max": DEFAULT_R0_T_MAX},
            verbose=verbose,
            prefix=prefix,
        )
        calc.variables.update(r0_fit_cache)
        r0_var = calc.get_variable(
            "r0",
            t_min=DEFAULT_R0_T_MIN,
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

    vprint("Calculating chi and F_chi...")
    all_chi = {}
    all_f_chi = {}
    r0_chi = None
    r0_chi_err = None
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

        for r_p, t_p in chi_pairs:
            try:
                chi_var = calc.get_variable("chi", R=r_p, T=t_p)
                chi_val = chi_var.get()
                if chi_val is not None and not np.isnan(chi_val):
                    all_chi[f"{r_p},{t_p}"] = float(chi_val)
            except Exception:
                continue

        r_for_force = sorted({p[0] for p in chi_pairs if p[1] == chi_t_large + 0.5})
        for r_p in r_for_force:
            try:
                f_var = calc.get_variable("F_chi", R=r_p, t_large=chi_t_large)
                f_val = f_var.get()
                if f_val is not None and not np.isnan(f_val):
                    all_f_chi[f"{r_p}"] = float(f_val)
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
        pass

    vprint("Calculating volume_r0...")
    volume_r0 = None
    volume_r0_err = None
    try:
        yaml_path = input_dir / "input.yaml"
        metro, _ = load_params(str(yaml_path))
        L_val = metro.L0

        if r0 is not None:
            vol_var = calc.get_variable(
                "volume_r0",
                L=L_val,
                t_min=DEFAULT_R0_T_MIN,
                t_max=DEFAULT_R0_T_MAX,
                target_force=sommer_target,
                r_min=DEFAULT_R0_R_MIN,
            )
            if vol_var.get() is not None and not np.isnan(vol_var.get()):
                volume_r0 = vol_var.get()
                volume_r0_err = vol_var.err()
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
        "a": float(lattice_spacing) if lattice_spacing is not None else None,
        "a_err": float(a_err) if a_err is not None else None,
        "tau_int": float(tau),
        "block_size": int(block_size),
        "analysis_settings": analysis_settings,
        "V_R": potentials,
        "V_R_err": potential_errors,
        "W_R_T": all_w,
        "r0_chi": float(r0_chi) if r0_chi is not None else None,
        "r0_chi_err": float(r0_chi_err) if r0_chi_err is not None else None,
        "chi": all_chi,
        "F_chi": all_f_chi,
        "creutz_P": all_creutz_P,
        "creutz_P_err": all_creutz_P_err,
        "a_creutz": all_a_creutz,
        "a_creutz_err": all_a_creutz_err,
        "plot_meta": {
            "potentials": potentials,
            "potential_errors": potential_errors,
            "cornell_params": cornell_params,
            "chi": all_chi,
            "F_chi": all_f_chi,
        },
    }


def get_or_calculate(
    path: str,
    force_recalc: bool = False,
    verbose: bool = True,
    calc_workers: Optional[int] = None,
    load_workers: int = 1,
    combine_equivalent_runs: bool = True,
) -> Dict[str, Any]:
    path = os.path.abspath(path)
    run_id = get_run_id(path)

    prefix = f"[{os.path.basename(path)}]"
    def vprint(msg: str):
        if verbose:
            import sys
            print(f"{prefix} {msg}", file=sys.stderr, flush=True)

    if not os.path.isdir(path):
        return {"error": "Directory not found"}

    if combine_equivalent_runs:
        analysis_id, grouped_paths = _discover_equivalent_runs(path)
    else:
        analysis_id = f"run_{run_id}"
        grouped_paths = [path]

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

        result = evaluate_run(
            combined_w_temp,
            Path(path),
            verbose=verbose,
            prefix=prefix,
            calc_workers=calc_workers,
        )
        result["path"] = path
        result["run_id"] = run_id
        result["analysis_id"] = analysis_id
        result["aggregation"] = aggregation

        if "error" not in result:
            save_result(analysis_id, result)
        return result
    except Exception as e:
        return {"error": str(e)}
