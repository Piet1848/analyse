import json
import hashlib
import os
from dataclasses import asdict
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple

import numpy as np

from load_input_yaml import load_params, MetropolisParams, GaugeObservableParams

import data_organizer as do
from calculator import Calculator

# --- CONFIGURATION ---
DATA_ROOT = Path("../data").resolve()
CALC_RESULT_BASE = DATA_ROOT / "calcResult"
THERMALIZATION_STEPS = 500
CALC_VERSION = "4.1"  # Bumped version: grouped run-combination before calculation
GROUP_IGNORE_METRO_FIELDS = {"seed", "nSweep"}

_GROUP_INDEX_CACHE: Dict[str, List[str]] | None = None
_GROUP_INDEX_ROOT: Path | None = None


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


def _load_combined_w_temp(run_paths: List[str]) -> Tuple[Optional[do.FileData], Dict[str, Any]]:
    w_temp_files: List[do.FileData] = []
    runs_with_w_temp = 0

    for run_path in run_paths:
        exp_data = do.ExperimentData(run_path)
        files = exp_data.data.get("W_temp", [])
        if files:
            runs_with_w_temp += 1
            w_temp_files.extend(files)

    combined = do.combine_file_data(
        w_temp_files,
        min_step=THERMALIZATION_STEPS,
        source_name="W_temp_combined",
    )

    combined_samples = 0
    if combined is not None and combined.observables:
        combined_samples = len(combined.observables[0].values)

    metadata = {
        "n_runs_in_group": len(run_paths),
        "n_runs_with_w_temp": runs_with_w_temp,
        "n_w_temp_files": len(w_temp_files),
        "n_samples_after_cut": combined_samples,
        "thermalization_steps": THERMALIZATION_STEPS,
    }
    return combined, metadata


def evaluate_run(file_data: do.FileData, input_dir: Path, sommer_target: float = 1.65) -> Dict[str, Any]:
    """
    Main analysis logic using Calculator.
    Expects preprocessed (thermalization-cut + combined) W_temp data.
    """
    if not file_data.observables:
        return {"error": "No combined W_temp observables available"}

    # 1. Calculate Autocorrelation (Tau) first
    calc_pre = Calculator(file_data)
    try:
        tau_var = calc_pre.get_variable("tau_int", obs_name="plaquette")
        tau = tau_var.get()
    except KeyError:
        tau = 0.5

    block_size = max(1, int(np.ceil(2 * tau)))

    # 2. Main Calculation (with Blocking)
    calc = Calculator(file_data, step_size=block_size)

    # --- A. Explicitly Calculate V(R) for all available R ---
    potentials = {}
    potential_errors = {}
    try:
        all_L = np.asarray(file_data.get("L").values)
        unique_Rs = sorted(np.unique(all_L))
        for r in unique_Rs:
            try:
                v_var = calc.get_variable("V_R", R=int(r))
                val = v_var.get()
                if val is not None and not np.isnan(val):
                    potentials[str(int(r))] = float(val)
                    potential_errors[str(int(r))] = float(v_var.err())
            except Exception:
                continue
    except Exception:
        pass

    # --- B. Collect W(R, T) for all available pairs ---
    all_w = {}
    try:
        all_T = np.array(file_data.get("T").values)
        pairs = sorted(list(set(zip(all_L, all_T))))
        for r_p, t_p in pairs:
            try:
                w_var = calc.get_variable("W_R_T", R=int(r_p), T=int(t_p))
                w_val = w_var.get()
                if w_val is not None:
                    all_w[f"{int(r_p)},{int(t_p)}"] = float(w_val)
            except Exception:
                continue
    except Exception:
        pass

    # --- C. Attempt r0 Calculation ---
    r0 = None
    r0_err = None
    cornell_params = None

    try:
        r0_var = calc.get_variable("r0", target_force=sommer_target)
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

    # --- D. Attempt Creutz Ratio (chi) analysis ---
    all_chi = {}
    all_f_chi = {}
    r0_chi = None
    r0_chi_err = None
    chi_t_large = 4

    try:
        chi_pairs = []
        if "all_L" in locals() and "all_T" in locals():
            r_values = sorted(np.unique(all_L))
            t_values = sorted(np.unique(all_T))
            for r in r_values:
                if r + 1 not in r_values:
                    continue
                for t in t_values:
                    if t + 1 not in t_values:
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

        r_for_force = sorted(list(set(p[0] for p in chi_pairs if p[1] == chi_t_large + 0.5)))
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
            discard_negative=True,
        )
        r0_chi = r0_chi_var.get()
        if np.isnan(r0_chi):
            r0_chi = None
        else:
            r0_chi_err = r0_chi_var.err()

    except Exception:
        pass

    # --- E. Volume Calculation ---
    volume_r0 = None
    volume_r0_err = None
    try:
        yaml_path = input_dir / "input.yaml"
        metro, _ = load_params(str(yaml_path))
        L_val = metro.L0

        if r0 is not None:
            vol_var = calc.get_variable("volume_r0", L=L_val, target_force=sommer_target)
            if vol_var.get() is not None and not np.isnan(vol_var.get()):
                volume_r0 = vol_var.get()
                volume_r0_err = vol_var.err()
    except Exception:
        pass

    # --- F. Creutz P-Ratio and a_creutz ---
    all_creutz_P = {}
    all_creutz_P_err = {}
    all_a_creutz = {}
    all_a_creutz_err = {}

    try:
        if "all_L" in locals():
            unique_Rs = sorted(np.unique(all_L))
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


def get_or_calculate(path: str, force_recalc: bool = False) -> Dict[str, Any]:
    path = os.path.abspath(path)
    run_id = get_run_id(path)

    if not os.path.isdir(path):
        return {"error": "Directory not found"}

    analysis_id, grouped_paths = _discover_equivalent_runs(path)

    if not force_recalc:
        cached = load_cached_result(analysis_id)
        if cached:
            result = dict(cached)
            result["path"] = path
            result["run_id"] = run_id
            result["analysis_id"] = analysis_id
            return result

    try:
        combined_w_temp, aggregation = _load_combined_w_temp(grouped_paths)
        if combined_w_temp is None:
            return {"error": "No W_temp data found in equivalent run group"}

        result = evaluate_run(combined_w_temp, Path(path))
        result["path"] = path
        result["run_id"] = run_id
        result["analysis_id"] = analysis_id
        result["aggregation"] = aggregation

        if "error" not in result:
            save_result(analysis_id, result)
        return result
    except Exception as e:
        return {"error": str(e)}
