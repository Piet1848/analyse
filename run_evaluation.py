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
CALC_VERSION = "2.0"  # Bump this to invalidate old caches automatically

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
                # Auto-invalidate if version mismatches
                if data.get("version") == CALC_VERSION:
                    return data
        except (json.JSONDecodeError, IOError): 
            return None
    return None

def save_result(run_id: str, data: Dict[str, Any]):
    CALC_RESULT_BASE.mkdir(parents=True, exist_ok=True)
    data["version"] = CALC_VERSION  # Stamp with version
    p = get_result_path(run_id)
    with open(p, 'w') as f: json.dump(data, f, indent=2)

def evaluate_run(data: do.ExperimentData, sommer_target: float = 1.65) -> Dict[str, Any]:
    """
    Main analysis logic using Calculator.
    """
    # 1. Get the main file (W_temp)
    if "W_temp" not in data.data or not data.data["W_temp"]:
        return {"error": "No W_temp data found"}
    
    file_data = data.data["W_temp"][0]
    file_data.remove_thermalization(THERMALIZATION_STEPS)
    
    if not file_data.observables:
        return {"error": "No data left after thermalization cut"}

    # 2. Calculate Autocorrelation (Tau) first
    # We use a default calculator (step=1) for this
    calc_pre = Calculator(file_data)
    try:
        # Try to find a good observable for autocorrelation
        tau_var = calc_pre.get_variable("tau_int", obs_name="plaquette")
        tau = tau_var.get()
    except KeyError:
        # Fallback to W_temp or just default
        tau = 0.5

    block_size = max(1, int(np.ceil(2 * tau)))

    # 3. Main Calculation (with Blocking)
    # Instantiate a NEW calculator with the correct block size
    calc = Calculator(file_data, step_size=block_size)
    try:
        # Calculate r0 (this triggers V_R calculations internally)
        r0_var = calc.get_variable("r0", target_force=sommer_target)
        r0 = r0_var.get()
        r0_err = r0_var.err()
    except Exception as e:
        return {"error": f"Calculation failed: {e}"}

    if r0 is None or np.isnan(r0):
        return {"error": "Could not determine Sommer parameter"}

    lattice_spacing = 0.5 / r0
    # Approx error propagation for a = 0.5/r0
    a_err = (0.5 / (r0**2)) * r0_err if r0_err else 0.0

    # 4. Harvest Plotting Metadata
    # We retrieve the V(R) values that were cached inside the calculator
    potentials = {}
    potential_errors = {}
    
    for key, var in calc.variables.items():
        if key[0] == "V_R":
            # key is ("V_R", frozenset({('R', 2), ...}))
            params = dict(key[1])
            R = params.get('R')
            if R is not None:
                potentials[R] = var.get()
                potential_errors[R] = var.err()
    
    # Retrieve fit parameters stored in Step 1
    cornell_params = r0_var.parameters.get("cornell_params", None)

    return {
        "r0": r0,
        "r0_err": r0_err,
        "a": lattice_spacing,
        "a_err": a_err,
        "tau_int": tau,
        "block_size": block_size,
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

        if "error" not in result:
            save_result(run_id, result)
        return result
    except Exception as e:
        return {"error": str(e)}