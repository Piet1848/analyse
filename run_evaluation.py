import json
import hashlib
import os
from pathlib import Path
from typing import Dict, Any, Optional
import numpy as np

# Import the new tools
import data_organizer as do
from calculator import Calculator

# --- CONFIGURATION ---
CALC_RESULT_BASE = Path("../data/calcResult")
THERMALIZATION_STEPS = 500 
CALC_VERSION = "2.1"  # Bumped version to invalidate old caches

def get_run_id(path: str) -> str:
    rel_path = os.path.relpath(path, start="../data") 
    return hashlib.md5(rel_path.encode('utf-8')).hexdigest()

def get_result_path(run_id: str) -> Path:
    return CALC_RESULT_BASE / f"{run_id}.json"

def load_cached_result(run_id: str) -> Optional[Dict[str, Any]]:
    p = get_result_path(run_id)
    if p.exists():
        try:
            with open(p, 'r') as f: 
                data = json.load(f)
                if data.get("version") == CALC_VERSION:
                    return data
        except (json.JSONDecodeError, IOError): 
            return None
    return None

def save_result(run_id: str, data: Dict[str, Any]):
    CALC_RESULT_BASE.mkdir(parents=True, exist_ok=True)
    data["version"] = CALC_VERSION
    p = get_result_path(run_id)
    with open(p, 'w') as f: json.dump(data, f, indent=2)

def evaluate_run(data: do.ExperimentData, sommer_target: float = 1.65) -> Dict[str, Any]:
    """
    Main analysis logic using Calculator. 
    Calculates as much as possible without failing completely if r0 is missing.
    """
    # 1. Get the main file (W_temp)
    if "W_temp" not in data.data or not data.data["W_temp"]:
        return {"error": "No W_temp data found"}
    
    file_data = data.data["W_temp"][0]
    file_data.remove_thermalization(THERMALIZATION_STEPS)
    
    if not file_data.observables:
        return {"error": "No data left after thermalization cut"}

    # 2. Calculate Autocorrelation (Tau) first
    calc_pre = Calculator(file_data)
    try:
        tau_var = calc_pre.get_variable("tau_int", obs_name="plaquette")
        tau = tau_var.get()
    except KeyError:
        tau = 0.5

    block_size = max(1, int(np.ceil(2 * tau)))

    # 3. Main Calculation (with Blocking)
    calc = Calculator(file_data, step_size=block_size)
    
    # --- A. Explicitly Calculate V(R) for all available R ---
    # We do this first so we have V_R even if r0 fails later.
    potentials = {}
    potential_errors = {}
    
    try:
        # Assuming L column corresponds to spatial separation R
        all_L = np.array(file_data.get("L").values)
        unique_Rs = sorted(np.unique(all_L))
        
        for r in unique_Rs:
            try:
                # Calculate V_R for this R
                v_var = calc.get_variable("V_R", R=int(r))
                val = v_var.get()
                if val is not None and not np.isnan(val):
                    potentials[int(r)] = val
                    potential_errors[int(r)] = v_var.err()
            except Exception:
                continue
    except Exception as e:
        # If we can't determine R values, we can't calculate V_R
        pass

    # --- B. Attempt r0 Calculation ---
    r0 = None
    r0_err = None
    cornell_params = None
    
    try:
        # This will use the cached V_R values from step A
        r0_var = calc.get_variable("r0", target_force=sommer_target)
        r0 = r0_var.get()
        # Treat NaN as None/Failure for r0
        if np.isnan(r0):
            r0 = None
        else:
            r0_err = r0_var.err()
            cornell_params = r0_var.parameters.get("cornell_params", None)
    except Exception:
        # r0 calculation failed (e.g., not enough points), but we proceed
        pass

    # Calculate lattice spacing if r0 was found
    lattice_spacing = None
    a_err = None
    if r0 is not None:
        lattice_spacing = 0.5 / r0
        a_err = (0.5 / (r0**2)) * r0_err if r0_err else 0.0

    return {
        "r0": r0,
        "r0_err": r0_err,
        "a": lattice_spacing,
        "a_err": a_err,
        "tau_int": tau,
        "block_size": block_size,
        "V_R": potentials, # Added top-level key for search_data.py
        "plot_meta": {
            "potentials": potentials,
            "potential_errors": potential_errors,
            "cornell_params": cornell_params
        }
    }

def get_or_calculate(path: str, force_recalc: bool = False) -> Dict[str, Any]:
    run_id = get_run_id(path)
    
    if not force_recalc:
        cached = load_cached_result(run_id)
        if cached: return cached

    try:
        if not os.path.isdir(path): return {"error": "Directory not found"}
        exp_data = do.ExperimentData(path)
        result = evaluate_run(exp_data)
        result["path"] = path
        result["run_id"] = run_id

        # Only save if we didn't hit a fatal error (missing files)
        # Partial results (r0=None) are fine to save/return
        if "error" not in result:
            save_result(run_id, result)
        return result
    except Exception as e:
        return {"error": str(e)}