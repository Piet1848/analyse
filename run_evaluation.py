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
CALC_VERSION = "2.6"  # Bumped version to invalidate old caches

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
    potentials = {}
    potential_errors = {}
    try:
        all_L = np.array(file_data.get("L").values)
        unique_Rs = sorted(np.unique(all_L))
        for r in unique_Rs:
            try:
                v_var = calc.get_variable("V_R", R=int(r))
                val = v_var.get()
                if val is not None and not np.isnan(val):
                    # Store as string keys and native floats
                    potentials[str(int(r))] = float(val)
                    potential_errors[str(int(r))] = float(v_var.err())
            except Exception: continue
    except Exception: pass

    # --- NEW: Collect W(R, T) for all available pairs ---
    all_w = {}
    try:
        all_T = np.array(file_data.get("T").values)
        # Find all unique (L, T) pairs
        pairs = sorted(list(set(zip(all_L, all_T))))
        for r_p, t_p in pairs:
            try:
                w_var = calc.get_variable("W_R_T", R=int(r_p), T=int(t_p))
                w_val = w_var.get()
                if w_val is not None:
                    all_w[f"{int(r_p)},{int(t_p)}"] = float(w_val)
            except: continue
    except: pass
    
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

    # --- C. Attempt Creutz Ratio (chi) analysis ---
    all_chi = {}
    all_f_chi = {}
    r0_chi = None
    r0_chi_err = None
    chi_t_large = 4 # Default large T for extracting the force

    try:
        # Determine R and T values for chi from Wilson loop pairs
        # chi is at R+0.5, T+0.5
        chi_pairs = []
        # Check if all_L and all_T were defined earlier
        if 'all_L' in locals() and 'all_T' in locals():
            r_values = sorted(np.unique(all_L))
            t_values = sorted(np.unique(all_T))
            for r in r_values:
                if r + 1 not in r_values: continue
                for t in t_values:
                    if t + 1 not in t_values: continue
                    chi_pairs.append((r + 0.5, t + 0.5))
        
        for r_p, t_p in chi_pairs:
            try:
                chi_var = calc.get_variable("chi", R=r_p, T=t_p)
                chi_val = chi_var.get()
                if chi_val is not None and not np.isnan(chi_val):
                    all_chi[f"{r_p},{t_p}"] = float(chi_val)
            except Exception:
                continue

        # Extract force F_chi at a fixed t_large
        r_for_force = sorted(list(set(p[0] for p in chi_pairs if p[1] == chi_t_large + 0.5)))
        for r_p in r_for_force:
             try:
                f_var = calc.get_variable("F_chi", R=r_p, t_large=chi_t_large)
                f_val = f_var.get()
                if f_val is not None and not np.isnan(f_val):
                    all_f_chi[f"{r_p}"] = float(f_val)
             except Exception:
                continue

        # Calculate r0 from chi
        r0_chi_var = calc.get_variable("r0_chi", t_large=chi_t_large, target_force=sommer_target, discard_negative=True)
        r0_chi = r0_chi_var.get()
        if np.isnan(r0_chi):
            r0_chi = None
        else:
            r0_chi_err = r0_chi_var.err()

    except Exception:
        # Chi analysis can fail, proceed without it
        pass
    
    # --- D. Creutz P-Ratio and a_creutz ---
    all_creutz_P = {}
    all_creutz_P_err = {}
    all_a_creutz = {}
    all_a_creutz_err = {}

    try:
        if 'all_L' in locals():
            unique_Rs = sorted(np.unique(all_L))
            even_Rs = [r for r in unique_Rs if r % 2 == 0 and r > 0] # Ensure R > 0

            for r in even_Rs:
                try:
                    # P-Ratio
                    p_var = calc.get_variable("creutz_P", R=int(r))
                    p_val = p_var.get()
                    if p_val is not None and not np.isnan(p_val):
                        all_creutz_P[str(int(r))] = float(p_val)
                        if p_var.err() is not None:
                            all_creutz_P_err[str(int(r))] = float(p_var.err())

                    # a_creutz
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