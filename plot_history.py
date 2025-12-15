#!/usr/bin/env python3
"""
Plot the thermalization/history of a specific observable.
Usage: python plot_history.py <path_to_run> <observable_name>
"""

import argparse
import sys
import matplotlib
# Try to use an interactive backend if available (specifically for WSL)
try:
    matplotlib.use('TkAgg')
except:
    pass

import matplotlib.pyplot as plt
import data_organizer as do

def main():
    parser = argparse.ArgumentParser(description="Plot observable history (thermalization).")
    parser.add_argument("path", help="Path to the run directory")
    parser.add_argument("observable", nargs="?", help="Name of the observable to plot")
    parser.add_argument("--xlim", type=int, help="Limit x-axis to first N steps", default=None)
    args = parser.parse_args()

    # 1. Load Data
    print(f"Loading data from {args.path}...")
    try:
        exp_data = do.ExperimentData(args.path)
    except Exception as e:
        sys.exit(f"Error loading data: {e}")

    if not exp_data.data:
        sys.exit(f"No .out files found in {args.path}")

    # 2. Inspect available observables from ALL files
    available_obs = set()
    for file_list in exp_data.data.values():
        for file_data in file_list:
            for obs in file_data.observables:
                available_obs.add(obs.name)
    
    sorted_obs = sorted(list(available_obs))

    # If user didn't provide an observable, list them and exit
    if not args.observable:
        print("\nPlease specify an observable to plot.")
        print("Available observables found in data:")
        for name in sorted_obs:
            print(f"  - {name}")
        sys.exit(1)

    target = args.observable
    if target not in available_obs:
        print(f"\nError: Observable '{target}' not found.")
        print(f"Available: {sorted_obs}")
        sys.exit(1)

    # 3. Plotting
    plt.figure(figsize=(12, 6))
    
    found_any = False
    
    # Iterate over all files to find the one containing the target
    for file_stem, file_list in exp_data.data.items():
        for file_data in file_list:
            # Check if this file contains the target observable
            target_obs = next((o for o in file_data.observables if o.name == target), None)
            
            # We also need a 'step' column in the SAME file to plot against
            step_obs = next((o for o in file_data.observables if o.name in ["# step", "step"]), None)

            if step_obs and target_obs:
                found_any = True
                
                # Align lengths (handle cases where a run crashed mid-write)
                min_len = min(len(step_obs.values), len(target_obs.values))
                steps = step_obs.values[:min_len]
                values = target_obs.values[:min_len]

                plt.plot(steps, values, '.', markersize=5, alpha=0.8, label=f"{file_stem}")

    if not found_any:
        sys.exit(f"Could not extract data for '{target}'. (It might be missing a corresponding 'step' column in its file).")

    plt.xlabel("Monte Carlo Step")
    plt.ylabel(target)
    plt.title(f"Thermalization History: {target}")
    plt.grid(True, linestyle='--', alpha=0.5)
    
    if args.xlim:
        plt.xlim(0, args.xlim)

    # Only show legend if we found multiple sources
    if len(exp_data.data) > 1:
        plt.legend()

    output_file = f"history_{target}.png"
    plt.savefig(output_file, dpi=150)
    print(f"Plot saved to {output_file}")
    
    # SHOW THE PLOT
    try:
        plt.show()
    except Exception as e:
        print(f"Warning: Could not display plot window ({e}). Plot is saved to file.")

if __name__ == "__main__":
    main()