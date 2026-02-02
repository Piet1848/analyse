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
from calculator import cornell_potential_ansatz
from load_input_yaml import load_params

# --- Configuration for Units ---
UNIT_MAP = {
    "a": "[fm]",
    "r0": "[fm]",
    "a_creutz": "[fm]",
    "creutz_P": "",
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
    # 1. Direct match in calc results (for static calculated fields)
    if calc and name in calc:
        return calc[name]

    # 2. Handle Dynamic tokens (e.g. V_R4, a_creutz2)
    # Check for error suffix first
    lookup_name = name
    is_err = False
    if name.endswith("_err"):
        lookup_name = name[:-4]
        is_err = True

    base, params = search_data.parse_dynamic_token(lookup_name)
    if base and calc:
        if base == "V_R":
            field = "V_R_err" if is_err else "V_R"
            return calc.get(field, {}).get(str(params["R"]))
        elif base == "creutz_P":
            field = "creutz_P_err" if is_err else "creutz_P"
            return calc.get(field, {}).get(str(params["R"]))
        elif base == "a_creutz":
            field = "a_creutz_err" if is_err else "a_creutz"
            return calc.get(field, {}).get(str(params["R"]))
        # Add other dynamic fields if needed (W_R_T etc)

    # 3. Metadata
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

def plot_potential_fit(matches: List[str]):
    """
    Specialized plotting function for V(r) vs r.
    Automatically detects varying parameters for the legend 
    and constant parameters for the title.
    Sorted by beta.
    """
    set_plot_style()
    fig, ax = plt.subplots()
    
    # --- 1. Pre-scan data to determine Legend vs Title ---
    valid_entries = []
    param_tracker = defaultdict(set)
    interesting_params = ["beta", "epsilon1", "epsilon2"]

    print(f"Loading metadata for {len(matches)} runs...", file=sys.stderr)

    for path in matches:
        try:
            yaml_path = os.path.join(path, "input.yaml")
            metro, gauge = load_params(yaml_path)
            
            calc_res = run_evaluation.get_or_calculate(path)
            
            # Recalc if plot_meta is missing
            if "error" not in calc_res and "plot_meta" not in calc_res:
                calc_res = run_evaluation.get_or_calculate(path, force_recalc=True)

            if "error" in calc_res or "plot_meta" not in calc_res:
                continue

            entry = {
                "path": path,
                "metro": metro,
                "meta": calc_res["plot_meta"]
            }
            valid_entries.append(entry)

            for p in interesting_params:
                val = getattr(metro, p)
                if isinstance(val, float):
                    val = round(val, 6)
                param_tracker[p].add(val)

        except Exception as e:
            print(f"Skipped {path}: {e}", file=sys.stderr)
            continue

    if not valid_entries:
        print("No valid potential data found to plot.", file=sys.stderr)
        return

    # --- NEW: Sort entries by beta ---
    # This ensures the legend appears as 2.2, 2.3, 2.4, ... instead of random order
    valid_entries.sort(key=lambda x: x["metro"].beta)

    # --- 2. Determine Labels ---
    legend_keys = [k for k in interesting_params if len(param_tracker[k]) > 1]
    title_keys = [k for k in interesting_params if len(param_tracker[k]) == 1]

    # --- 3. Plotting Loop ---
    for entry in valid_entries:
        meta = entry["meta"]
        metro = entry["metro"]
        
        raw_pots = meta["potentials"]
        raw_errs = meta["potential_errors"]
        pots = {int(k): v for k, v in raw_pots.items()}
        errs = {int(k): v for k, v in raw_errs.items()}
        params = meta["cornell_params"]
        
        if not pots: 
            continue

        rs = sorted(pots.keys())
        Vs = [pots[r] for r in rs]
        V_errs = [errs.get(r, 0.0) for r in rs]
        
        # Generate Dynamic Label
        if legend_keys:
            label_parts = [f"{k}={getattr(metro, k)}" for k in legend_keys]
            label = ", ".join(label_parts)
        else:
            label = os.path.basename(entry["path"])

        # Plot Data
        p = ax.errorbar(rs, Vs, yerr=V_errs, fmt='o', label=label, capsize=3)
        color = p[0].get_color()
        
        # Plot Fit
        if params:
            r_dense = np.linspace(min(rs), max(rs), 100)
            V_fit = cornell_potential_ansatz(
               r_dense, params['A'], params['sigma'], params['B']
            )
            ax.plot(r_dense, V_fit, '--', color=color, alpha=0.7)

    # --- 4. Titles and Formatting ---
    base_title = "Static Quark Potential V(r)"
    if title_keys:
        subset_str = ", ".join([f"{k}={list(param_tracker[k])[0]}" for k in title_keys])
        ax.set_title(f"{base_title}\n({subset_str})")
    else:
        ax.set_title(base_title)

    ax.set_xlabel("r/a")
    ax.set_ylabel("aV(r)")
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.3)
    
    plot_file = "plot_potential.png"
    plt.tight_layout()
    plt.savefig(plot_file, dpi=150)
    print(f"Plot saved to {plot_file}")
    plt.show()

def plot_creutz_data(matches: List[str]):
    """
    Specialized plotting function for Creutz Ratio analysis.
    - Plots chi(R,T) vs T for various R.
    - Plots R^2 * F(R) vs R to show r0 determination.
    """
    set_plot_style()
    # Using 2 subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 12), sharex=False)
    
    # --- 1. Pre-scan and data collection (similar to plot_potential_fit) ---
    valid_entries = []
    param_tracker = defaultdict(set)
    interesting_params = ["beta", "epsilon1", "epsilon2"]

    print(f"Loading metadata for {len(matches)} runs...", file=sys.stderr)

    for path in matches:
        try:
            yaml_path = os.path.join(path, "input.yaml")
            metro, gauge = load_params(yaml_path)
            
            calc_res = run_evaluation.get_or_calculate(path, force_recalc=False)
            if "error" in calc_res or "plot_meta" not in calc_res or "chi" not in calc_res["plot_meta"]:
                continue

            entry = {"path": path, "metro": metro, "meta": calc_res["plot_meta"]}
            valid_entries.append(entry)

            for p in interesting_params:
                param_tracker[p].add(getattr(metro, p, None))

        except Exception as e:
            print(f"Skipped {path}: {e}", file=sys.stderr)
            continue
    
    if not valid_entries:
        print("No valid Creutz ratio data found to plot.", file=sys.stderr)
        return

    valid_entries.sort(key=lambda x: x["metro"].beta)
    legend_keys = [k for k in interesting_params if len(param_tracker[k]) > 1]
    title_keys = [k for k in interesting_params if len(param_tracker[k]) == 1]

    # --- 2. Plotting Loop ---
    markers = ['o', 's', '^', 'D', 'v', '<', '>']
    
    for i, entry in enumerate(valid_entries):
        meta = entry["meta"]
        metro = entry["metro"]
        
        all_chi = meta.get("chi", {})
        f_chi = meta.get("F_chi", {})

        if not all_chi:
            continue

        # Generate Dynamic Label
        if legend_keys:
            label = ", ".join([f"{k}={getattr(metro, k)}" for k in legend_keys])
        else:
            label = os.path.basename(entry["path"])
        
        color = plt.cm.viridis(i / max(1, len(valid_entries)))
        
        # --- Plot 1: chi(R,T) vs T ---
        # Group chi values by R
        chi_by_r = defaultdict(list)
        for key, val in all_chi.items():
            r, t = map(float, key.split(','))
            chi_by_r[r].append((t, val))
        
        # Plot for a few R values to avoid clutter
        sorted_rs = sorted(chi_by_r.keys())
        for r_idx, r_val in enumerate(sorted_rs[:4]): # Plot for first 4 R values
            ts, chis = zip(*sorted(chi_by_r[r_val]))
            marker = markers[r_idx % len(markers)]
            ax1.plot(ts, chis, marker=marker, linestyle='-', color=color, label=f'{label}, R={r_val}')

        # --- Plot 2: R^2 * F(R) vs R ---
        if f_chi:
            rs = sorted([float(k) for k in f_chi.keys()])
            r2_f = [r**2 * f_chi[str(r)] for r in rs]
            ax2.plot(rs, r2_f, marker='o', linestyle='--', color=color, label=label)

    # --- 3. Final Touches ---
    # Plot 1 Formatting
    ax1.set_xlabel("T/a")
    ax1.set_ylabel("a$^2 \chi(R,T)$")
    ax1.set_title("Creutz Ratios vs. Time")
    ax1.legend(fontsize='small', ncol=max(1, len(valid_entries)))
    ax1.grid(True, linestyle='--', alpha=0.3)

    # Plot 2 Formatting
    ax2.axhline(1.65, color='red', linestyle=':', label="Sommer r0 Target (1.65)")
    ax2.set_xlabel("R/a")
    ax2.set_ylabel("R$^2$ a$^2$F(R)")
    ax2.set_title("Static Force from Creutz Ratios")
    ax2.legend(fontsize='small')
    ax2.grid(True, linestyle='--', alpha=0.3)
    
    # Global Title
    base_title = "Creutz Ratio Analysis"
    if title_keys:
        subset_str = ", ".join([f"{k}={list(param_tracker[k])[0]}" for k in title_keys])
        fig.suptitle(f"{base_title}\n({subset_str})", fontsize=18)
    else:
        fig.suptitle(base_title, fontsize=18)

    plot_file = "plot_creutz.png"
    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) # Adjust for suptitle
    plt.savefig(plot_file, dpi=150)
    print(f"Plot saved to {plot_file}")
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Search and Plot Lattice Data")
    parser.add_argument("root", help="Root data directory")
    # Make tokens optional so we can just run a search without x/y if needed
    parser.add_argument("tokens", nargs="*", help="Filters (param=val) and Plot Config (x=..., y=...)")
    parser.add_argument("--no-agg", action="store_true", help="Disable aggregation of duplicate points")
    # NEW ARGUMENT
    parser.add_argument("--plot-potential", action="store_true", help="Plot V(r) fit instead of scalar trends")
    parser.add_argument("--plot-creutz", action="store_true", help="Plot Creutz ratio analysis (chi vs T, F vs R)")
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
            # Allow clean pass-through if user didn't specify x/y but we don't need them for potential
            pass

    # Validation: If NOT plotting potential or creutz, we strictly need X and Y
    if not args.plot_potential and not args.plot_creutz:
        if not plot_config["x"] or not plot_config["y"]:
            print("Error: You must specify x and y axes (e.g., x=beta y=r0) OR use --plot-potential/--plot-creutz", file=sys.stderr)
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

    print(f"Found {len(matches)} matching runs.", file=sys.stderr)

    # --- BRANCH: Plot Potential vs Plot Trend ---
    if args.plot_potential:
        plot_potential_fit(matches)
        sys.exit(0)
    elif args.plot_creutz:
        plot_creutz_data(matches)
        sys.exit(0)

    # --- 3. Collect Data (Standard Mode) ---
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

    # --- 4. Plotting (Standard Mode) ---
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