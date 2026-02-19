import numpy as np
import matplotlib.pyplot as plt

# Local modules
import data_organizer as do
from calculator import Calculator, cornell_potential_ansatz

# Setup plotting style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['font.size'] = 12

# --- CONFIGURATION ---

# 1. Data Path 
# Example: "../data/20251221/45" 
FILE_PATH = "../data/20260205/36"

# 2. Analysis Parameters
THERMALIZATION_STEPS = 500
SOMMER_TARGET_FORCE = 1.65

print(f"Analyzing Directory: {FILE_PATH}")
print(f"Parameters: Therm={THERMALIZATION_STEPS}, Target={SOMMER_TARGET_FORCE}")

# --- 1. Load Data and Initialize Calculator ---

try:
    # Initialize Experiment Data
    experiment = do.ExperimentData(FILE_PATH)
    
    if "W_temp" not in experiment.data or not experiment.data["W_temp"]:
        print("Error: No 'W_temp.out' file found in the directory.")
    else:
        # Get the first W_temp file (assuming only one usually or we take the first)
        w_temp_data = experiment.data["W_temp"][0]
        print(f"Loaded W_temp data with {len(w_temp_data.observables[0].values)} entries (after thermalization cut if applied).")
        calc = Calculator(w_temp_data, n_bootstrap=50, step_size=1)
        print("Calculator initialized.")

        # --- 2. Calculate and Display Results ---
        r0_var = calc.get_variable("r0", t_min=5, target_force=1.65)
        print(f"Calculated r0/a: {r0_var.value:.3f} ± {r0_var.error:.3f}")

        # --- 3. Plotting ---
        try:
            # Get V(R) data points
            unique_L = sorted(np.unique(calc.file_data.get("L").values))
            rs = []
            vs = []
            vs_err = []

            for r in unique_L:
                try:
                    # Retrieve the cached V(R) variable used for r0 calculation
                    v_r = calc.get_variable("V_R", R=int(r), t_min=5)
                    if v_r.value is not None and not np.isnan(v_r.value):
                        rs.append(r)
                        vs.append(v_r.value)
                        vs_err.append(v_r.error if v_r.error is not None else 0)
                except Exception:
                    continue

            rs = np.array(rs)
            vs = np.array(vs)
            vs_err = np.array(vs_err)

            plt.figure()
            # Plot measured V(R) points
            plt.errorbar(rs, vs, yerr=vs_err, fmt='o', label='Measured V(R)', capsize=3)

            # Plot Cornell Potential Fit
            if "cornell_params" in r0_var.parameters:
                params = r0_var.parameters["cornell_params"]
                A, sigma, B = params['A'], params['sigma'], params['B']

                # Generate smooth curve for the fit
                x_fit = np.linspace(min(rs), max(rs), 100)
                y_fit = cornell_potential_ansatz(x_fit, A, sigma, B)
                
                plt.plot(x_fit, y_fit, '-', label=f'Fit: $\sigma={sigma:.3f}, B={B:.3f}$')
                
                # Mark r0 on the plot
                r0_val = r0_var.value
                if r0_val is not None:
                    # V(r0) point
                    v_r0 = cornell_potential_ansatz(r0_val, A, sigma, B)
                    plt.plot(r0_val, v_r0, 'r*', markersize=12, label=f'$r_0/a = {r0_val:.3f}$')
                    plt.axvline(x=r0_val, color='r', linestyle='--', alpha=0.5)

            plt.xlabel('R/a')
            plt.ylabel('aV(R)')
            plt.title(f'Static Potential V(R) vs R (t_min=5)')
            plt.legend()
            plt.grid(True)
            plt.tight_layout()
            plt.show()
            print("Plot displayed.")

        except Exception as plot_err:
            print(f"An error occurred during plotting: {plot_err}")

except Exception as e:
    print(f"An error occurred loading data: {e}")