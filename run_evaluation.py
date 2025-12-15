import json
import hashlib
import os
from pathlib import Path
from typing import Dict, Any, Optional

import numpy as np
import data_organizer as do
import analyze_wilson as wilson
import argparse

CALC_RESULT_BASE = Path("../data/calcResult")

# --- CONFIGURATION ---
# Number of update steps to discard as thermalization
THERMALIZATION_STEPS = 500 
# ---------------------

def get_run_id(path: str) -> str:
    """
    Generate a unique, filesystem-safe ID for a dataset path.
    We use a hash of the absolute path to ensure uniqueness.
    """
    abs_path = os.path.abspath(path)
    return hashlib.md5(abs_path.encode('utf-8')).hexdigest()

def get_result_path(run_id: str) -> Path:
    return CALC_RESULT_BASE / f"{run_id}.json"

def load_cached_result(run_id: str) -> Optional[Dict[str, Any]]:
    """Try to load results from JSON."""
    p = get_result_path(run_id)
    if p.exists():
        try:
            with open(p, 'r') as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError):
            return None
    return None

def save_result(run_id: str, data: Dict[str, Any]):
    """Save results to JSON."""
    CALC_RESULT_BASE.mkdir(parents=True, exist_ok=True)
    p = get_result_path(run_id)
    with open(p, 'w') as f:
        json.dump(data, f, indent=2)

def sommer_parameter(data: do.ExperimentData, sommer_target: float = 1.65) -> Dict[str, Any]:
    """
    Calculate Sommer parameter and Lattice spacing.
    Uses iterative bootstrapping to save RAM.
    """
    if "W_temp" not in data.data or not data.data["W_temp"]:
        return {"error": "No W_temp data found"}

    fileData = data.data["W_temp"][0] 
    
    # 1. Apply Thermalization Cut
    fileData.remove_thermalization(THERMALIZATION_STEPS)
    
    if not fileData.observables or len(fileData.observables[0].values) == 0:
        return {"error": "No data left after thermalization cut"}

    sommer_parameters = []
    lattice_spacings = []
    
    # --- SMART MEMORY FIX ---
    # Instead of generating 200 samples at once, we loop and generate 1 at a time.
    # This keeps RAM usage constant (low) throughout the process.
    n_bootstrap_samples = 200
    base_seed = 42

    for i in range(n_bootstrap_samples):
        # Generate exactly ONE bootstrap sample
        # We use (base_seed + i) to ensure every sample is unique but reproducible
        # [0] is needed because get_bootstrap returns a list
        fd = fileData.get_bootstrap(n_bootstrap=1, seed=base_seed + i)[0]

        # --- Analysis Logic (Same as before) ---
        records = wilson.load_w_temp(fd.observables)
        averages = wilson.average_wilson_loops(records)
        
        # Fit Potential
        potentials, errors = wilson.fit_potential_from_time(averages, t_min=2)
        
        # Fit Sommer Parameter
        r0_over_a, _ = wilson.fit_sommer_parameter(
            potentials, 
            errors=errors,
            target_force_r2=sommer_target
        )

        if r0_over_a and r0_over_a > 0:
            lattice_spacing = 0.5 / r0_over_a 
            sommer_parameters.append(r0_over_a)
            lattice_spacings.append(lattice_spacing)
        
        # Python automatically frees the memory for 'fd' here because it goes out of scope
    
    if not lattice_spacings:
        return {"error": "Could not determine Sommer parameter"}

    mean_a = float(np.mean(lattice_spacings))
    stddev_a = float(np.std(lattice_spacings, ddof=1))
    mean_r0 = float(np.mean(sommer_parameters))
    stddev_r0 = float(np.std(sommer_parameters, ddof=1))

    return {
        "r0": mean_r0,
        "r0_err": stddev_r0,
        "a": mean_a,
        "a_err": stddev_a,
        "n_samples": len(lattice_spacings),
        # --- NEW: Return nominal data for plotting ---
        "plot_meta": {
            "potentials": pots_nominal,          # The V(r) points
            "potential_errors": errs_nominal,    # Errors on V(r)
            "cornell_params": params_nominal     # The fit A, B, sigma
        }
    }

def get_or_calculate(path: str, force_recalc: bool = False) -> Dict[str, Any]:
    """
    The main entry point for other tools.
    Checks cache -> loads data -> calculates -> saves -> returns.
    """
    run_id = get_run_id(path)
    
    if not force_recalc:
        cached = load_cached_result(run_id)
        if cached:
            return cached

    try:
        if not os.path.isdir(path):
             return {"error": "Directory not found"}

        exp_data = do.ExperimentData(path)

        result = sommer_parameter(exp_data)
        
        result["path"] = path
        result["run_id"] = run_id

        if "error" not in result:
            save_result(run_id, result)
            
        return result
        
    except Exception as e:
        return {"error": str(e)}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=str, help="Path to dataset folder")
    parser.add_argument("--force", action="store_true", help="Force recalculation")
    args = parser.parse_args()
    
    print(f"Processing {args.path}...")
    res = get_or_calculate(args.path, force_recalc=args.force)
    print(json.dumps(res, indent=2))