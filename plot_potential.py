#!/usr/bin/env python3
"""
Plot V(R) with weighted fitting and dynamic range adjustment.
"""

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt

import data_organizer as do
import analyze_wilson as wilson

def iterative_cornell_fit(potentials, errors, limit_r0_factor=None):
    """
    1. Weighted fit on ALL data.
    2. If limit_r0_factor is set, filter data where R > limit_r0_factor * r0.
    3. Refit on subset.
    """
    # 1. First Pass: Weighted Fit with all data
    r0_initial, params = wilson.fit_sommer_parameter(potentials, errors=errors)
    
    if r0_initial is None or not limit_r0_factor:
        return r0_initial, params

    # 2. Dynamic Filtering
    max_valid_R = r0_initial * limit_r0_factor
    
    print(f"  > Initial r0/a = {r0_initial:.3f}")
    print(f"  > Dynamic Cutoff: Excluding R > {max_valid_R:.2f} ...")

    filtered_potentials = {r: v for r, v in potentials.items() if r <= max_valid_R}
    
    # If we filtered too much, revert
    if len(filtered_potentials) < 3:
        print("  > Warning: Cutoff left too few points. Reverting to full set.")
        return r0_initial, params

    # 3. Second Pass: Refit
    r0_final, params_final = wilson.fit_sommer_parameter(filtered_potentials, errors=errors)
    
    if r0_final:
        print(f"  > Final r0/a = {r0_final:.3f}")
    else:
        print("  > Final fit failed, reverting.")
        return r0_initial, params
        
    return r0_final, params_final

def main():
    parser = argparse.ArgumentParser(description="Plot Static Potential V(R)")
    parser.add_argument("path", help="Path to the specific run directory")
    parser.add_argument("--tmin", type=int, default=2, help="Minimum t for potential fitting")
    # New argument to control dynamic range
    parser.add_argument("--cutoff", type=float, default=0.0, 
                        help="If > 0, exclude points where R > cutoff * r0 (e.g. 2.5)")
    args = parser.parse_args()

    # --- Load Data ---
    print(f"Loading data from {args.path}...")
    try:
        exp_data = do.ExperimentData(args.path)
    except Exception as e:
        sys.exit(f"Error loading data: {e}")

    if "W_temp" not in exp_data.data or not exp_data.data["W_temp"]:
        sys.exit("Error: No W_temp data found.")

    file_data = exp_data.data["W_temp"][0]
    records = wilson.load_w_temp(file_data.observables)
    averages = wilson.average_wilson_loops(records)
    
    # Get Potentials and Errors
    potentials, errors = wilson.fit_potential_from_time(averages, t_min=args.tmin)
    
    if not potentials:
        sys.exit("Could not fit any potentials.")

    # --- Clean Errors for Fitting ---
    # Convert errors that are None or 0 to a small epsilon or average
    # so they don't break the weighted fit
    cleaned_errors = {}
    valid_errs = [e for e in errors.values() if e > 0]
    mean_err = np.mean(valid_errs) if valid_errs else 0.01

    for r, err in errors.items():
        if err <= 0 or err is None:
            cleaned_errors[r] = mean_err * 100 # Give it tiny weight
        else:
            cleaned_errors[r] = err

    # --- Fit ---
    limit_factor = args.cutoff if args.cutoff > 0 else None
    r0, cornell_params = iterative_cornell_fit(potentials, cleaned_errors, limit_factor)

    # --- Plotting ---
    plt.figure(figsize=(10, 6))
    
    Rs = sorted(potentials.keys())
    Vs = [potentials[r] for r in Rs]
    Es = [cleaned_errors[r] for r in Rs]

    # Plot ALL data points (black)
    plt.errorbar(Rs, Vs, yerr=Es, fmt='o', capsize=3, color='black', label='Lattice Data $V(R)$')

    # Plot Fit
    if cornell_params:
        sigma = cornell_params['sigma']
        B = cornell_params['B']
        A = cornell_params['A']
        
        # Determine plot range
        max_R_plot = max(Rs)
        if r0 and r0 > max_R_plot: max_R_plot = r0 * 1.1
            
        r_smooth = np.linspace(min(Rs), max_R_plot, 200)
        v_smooth = wilson.cornell_potential_ansatz(r_smooth, A, sigma, B)
        
        plt.plot(r_smooth, v_smooth, '-', color='red', alpha=0.8, 
                 label=f'Cornell Fit\n$\\sigma={sigma:.3f}, B={B:.3f}$')
        
        if r0:
             plt.axvline(r0, color='blue', linestyle='--', alpha=0.5, label=f'$r_0/a = {r0:.3f}$')

        # Visual indicator of cutoff if used
        if limit_factor and r0:
            cutoff_val = r0 * limit_factor
            if cutoff_val < max(Rs):
                plt.axvspan(cutoff_val, max(Rs), color='gray', alpha=0.1, label='Excluded from Fit')

    plt.xlabel("Spatial Separation $R/a$")
    plt.ylabel("Static Potential $aV(R)$")
    plt.title(f"Static Potential for {args.path}")
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.legend()
    
    output_file = "potential_plot_weighted.png"
    plt.savefig(output_file, dpi=150)
    print(f"Plot saved to {output_file}")
    plt.show()

if __name__ == "__main__":
    main()