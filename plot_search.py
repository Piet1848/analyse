#!/usr/bin/env python3
"""
Search for data and plot the results immediately.
"""

import argparse
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from typing import Dict, Any, List, Tuple

# Import your existing tools
import search_data
import run_evaluation
from load_input_yaml import load_params

# --- Configuration for Units ---
UNIT_MAP = {
    "a": "[fm]",
    "r0": "[fm]",
    "r0_over_a": "",
    "beta": "",
    "sigma": "[fm$^{-2}$]"
}

# --- Styling Configuration ---
def set_plot_style():
    plt.rcParams.update({
        'font.size': 14,
        'axes.labelsize': 16,
        'axes.titlesize': 16,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'legend.fontsize': 14,
        'figure.figsize': (11, 7),
        'lines.linewidth': 1.2,
        'lines.markersize': 6,
        'errorbar.capsize': 3
    })

def get_value(name: str, metro: Any, gauge: Any, calc: Dict[str, Any]) -> Any:
    """Helper to extract a value from metadata or calculated results."""
    if calc and name in calc:
        return calc[name]
    if hasattr(metro, name):
        return getattr(metro, name)
    if hasattr(gauge, name):
        return getattr(gauge, name)
    return None

def weighted_mean(values: List[float], errors: List[float]) -> Tuple[float, float]:
    """Calculate weighted mean and propagated error."""
    vals = np.array(values)
    errs = np.array([e if e else 0.0 for e in errors])
    
    if np.any(errs <= 0):
        mean = np.mean(vals)
        err = np.std(vals, ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0.0
        return mean, err

    weights = 1.0 / (errs**2)
    w_sum = np.sum(weights)
    mean = np.sum(vals * weights) / w_sum
    err = 1.0 / np.sqrt(w_sum)
    return mean, err

def main():
    parser = argparse.ArgumentParser(description="Search and Plot Lattice Data")
    parser.add_argument("root", help="Root data directory")
    parser.add_argument("tokens", nargs="*", help="Filters (param=val) and Plot Config (x=..., y=..., hue=...)")
    parser.add_argument("--no-agg", action="store_true", help="Disable aggregation of duplicate points")
    args = parser.parse_args()

    # --- 1. Parse Arguments ---
    search_tokens = []
    plot_config = {"x": None, "y": None, "hue": None}
    
    for token in args.tokens:
        if "=" in token:
            key, val = token.split("=", 1)
            key = key.strip()
            if key in ["x", "y", "hue"]:
                plot_config[key] = val.strip()
            else:
                search_tokens.append(token)
        else:
            print(f"Ignored invalid token: {token}", file=sys.stderr)

    if not plot_config["x"] or not plot_config["y"]:
        print("Error: You must specify x and y axes (e.g., x=beta y=r0)", file=sys.stderr)
        sys.exit(1)

    # --- 2. Find Matching Runs ---
    print(f"Searching in {args.root}...", file=sys.stderr)
    try:
        criteria, _ = search_data.parse_tokens(search_tokens)
        matches = search_data.find_matching_runs(args.root, criteria)
    except Exception as e:
        print(f"Search Error: {e}", file=sys.stderr)
        sys.exit(1)

    if not matches:
        print("No matching runs found.", file=sys.stderr)
        sys.exit(0)

    print(f"Found {len(matches)} matching runs. Processing...", file=sys.stderr)

    # --- 3. Collect Data ---
    raw_data = defaultdict(list)
    x_name = plot_config["x"]
    y_name = plot_config["y"]
    hue_name = plot_config["hue"]
    y_err_name = f"{y_name}_err"

    for path in matches:
        try:
            metro, gauge = load_params(os.path.join(path, "input.yaml"))
            calc_res = run_evaluation.get_or_calculate(path)
            
            if "error" in calc_res:
                continue

            x_val = get_value(x_name, metro, gauge, calc_res)
            y_val = get_value(y_name, metro, gauge, calc_res)
            y_err = get_value(y_err_name, metro, gauge, calc_res)
            
            hue_val = "All"
            if hue_name:
                h = get_value(hue_name, metro, gauge, calc_res)
                if isinstance(h, float):
                    h = round(h, 6)
                hue_val = h if h is not None else "None"

            if x_val is not None and y_val is not None:
                raw_data[hue_val].append((x_val, y_val, y_err))
                
        except Exception as e:
            print(f"Error processing {path}: {e}", file=sys.stderr)
            continue

    # --- 4. Plotting ---
    set_plot_style()
    fig, ax = plt.subplots()
    
    # Sort hues
    sorted_hues = sorted(raw_data.keys(), key=lambda k: (str(k) if not isinstance(k, (int, float)) else k))
    
    markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*']

    for i, hue_val in enumerate(sorted_hues):
        points = raw_data[hue_val]
        
        # duplicates (same X)
        if not args.no_agg:
            points_by_x = defaultdict(list)
            for x, y, err in points:
                points_by_x[x].append((y, err))
            
            agg_points = []
            for x in sorted(points_by_x.keys()):
                ys = [p[0] for p in points_by_x[x]]
                es = [p[1] for p in points_by_x[x]]
                
                if len(ys) > 1:
                    mean, err = weighted_mean(ys, es)
                    agg_points.append((x, mean, err))
                else:
                    agg_points.append((x, ys[0], es[0]))
            points = agg_points
        else:
            points.sort(key=lambda p: p[0])

        xs = [p[0] for p in points]
        ys = [p[1] for p in points]
        yerrs = [p[2] for p in points]
        
        label = f"{hue_name}={hue_val}" if hue_name else y_name
        
        clean_errs = [e if e else 0 for e in yerrs]
        
        marker = markers[i % len(markers)]

        ax.errorbar(
            xs, ys, 
            yerr=clean_errs, 
            fmt=f'-{marker}',
            label=label, 
            alpha=0.9, 
            elinewidth=1.2,
            capthick=1.2
        )

    # --- Labels & Units ---
    x_label = f"{x_name} {UNIT_MAP.get(x_name, '')}"
    y_label = f"{y_name} {UNIT_MAP.get(y_name, '')}"

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    
    title = f"{y_name} vs {x_name}"
    if hue_name:
        title += f" (grouped by {hue_name})"
    ax.set_title(title)
    
    ax.grid(True, linestyle='--', alpha=0.3)
    
    if hue_name:
        ax.legend(title=hue_name, frameon=True, framealpha=0.9)
    
    plt.tight_layout()
    
    plot_file = "plot_output.png"
    plt.savefig(plot_file, dpi=150)
    print(f"Plot saved to {plot_file}")
    plt.show()

if __name__ == "__main__":
    main()