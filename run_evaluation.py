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
THERMALIZATION_STEPS = 500 

def get_run_id(path: str) -> str:
    abs_path = os.path.abspath(path)
    return hashlib.md5(abs_path.encode('utf-8')).hexdigest()

def get_result_path(run_id: str) -> Path:
    return CALC_RESULT_BASE / f"{run_id}.json"

def load_cached_result(run_id: str) -> Optional[Dict[str, Any]]:
    p = get_result_path(run_id)
    if p.exists():
        try:
            with open(p, 'r') as f: return json.load(f)
        except (json.JSONDecodeError, IOError): return None
    return None

def save_result(run_id: str, data: Dict[str, Any]):
    CALC_RESULT_BASE.mkdir(parents=True, exist_ok=True)
    p = get_result_path(run_id)
    with open(p, 'w') as f: json.dump(data, f, indent=2)

def sommer_parameter(data: do.ExperimentData, sommer_target: float = 1.65) -> Dict[str, Any]:
    # 1. Calculate Autocorrelation (tau_int)
    # We look for Plaquette or Retrace data to determine the global autocorrelation
    tau_int_val = 0.5
    
    obs_file = None
    target_obs_name = None
    
    if "plaquette" in data.data and data.data["plaquette"]:
        obs_file = data.data["plaquette"][0]
        target_obs_name = "plaquette"
    elif "retrace" in data.data and data.data["retrace"]:
        obs_file = data.data["retrace"][0]
        target_obs_name = "retrace"
        
    if obs_file:
        # Create a temp copy to cut thermalization without affecting original structure yet
        # (Though we are only reading values, so it's fine)
        # Find observables
        step_obs = next((o for o in obs_file.observables if o.name in ["# step", "step"]), None)
        val_obs = next((o for o in obs_file.observables if o.name == target_obs_name), None)
        
        if step_obs and val_obs:
            steps = np.array(step_obs.values)
            vals = np.array(val_obs.values)
            
            # Apply thermalization cut for calculation
            mask = steps >= THERMALIZATION_STEPS
            if np.sum(mask) > 100:
                clean_vals = vals[mask]
                tau_int_val = wilson.calculate_tau_int(clean_vals)

    # 2. Setup Bootstrapping with Blocking
    # Block size should be > 2 * tau_int
    block_size_steps = max(1, int(np.ceil(2 * tau_int_val)))

    if "W_temp" not in data.data or not data.data["W_temp"]:
        return {"error": "No W_temp data found"}

    fileData = data.data["W_temp"][0] 
    fileData.remove_thermalization(THERMALIZATION_STEPS)
    
    if not fileData.observables or len(fileData.observables[0].values) == 0:
        return {"error": "No data left after thermalization cut"}

    sommer_parameters = []
    lattice_spacings = []
    
    n_bootstrap_samples = 200
    base_seed = 42

    for i in range(n_bootstrap_samples):
        # Use the new Blocked Bootstrap
        # This samples STEPS, not rows, maintaining L/T correlations
        fd = fileData.get_blocked_bootstrap(
            n_bootstrap=1, 
            block_size=block_size_steps, 
            seed=base_seed + i
        )[0]

        records = wilson.load_w_temp(fd.observables)
        averages = wilson.average_wilson_loops(records)
        potentials, errors = wilson.fit_potential_from_time(averages, t_min=2)
        r0_over_a, _ = wilson.fit_sommer_parameter(potentials, errors=errors, target_force_r2=sommer_target)

        if r0_over_a and r0_over_a > 0:
            lattice_spacing = 0.5 / r0_over_a 
            sommer_parameters.append(r0_over_a)
            lattice_spacings.append(lattice_spacing)
    
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
         "tau_int": float(tau_int_val),
         "block_size": block_size_steps,
        "n_samples": len(lattice_spacings),
    }


def get_or_calculate(path: str, force_recalc: bool = False) -> Dict[str, Any]:
    run_id = get_run_id(path)
    
    if not force_recalc:
        cached = load_cached_result(run_id)
        if cached: return cached

    try:
        if not os.path.isdir(path): return {"error": "Directory not found"}
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