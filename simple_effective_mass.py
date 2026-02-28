import numpy as np
import data_organizer as do
from calculator import Calculator

def main():
    # Configuration
    # Adjust this path to point to your data directory relative to this script
    file_path = "../data/20260205/36" 
    R_target = 6
    
    print(f"Analyzing: {file_path}")
    
    try:
        # 1. Load Data
        # ExperimentData loads all observable files in the directory
        experiment = do.ExperimentData(file_path)
        
        if "W_temp" not in experiment.data or not experiment.data["W_temp"]:
            print("Error: W_temp data not found in the specified directory.")
            return
            
        # We use the first W_temp dataset found
        w_data = experiment.data["W_temp"][0]
        
        # 2. Initialize Calculator
        # n_bootstrap determines how many bootstrap samples are generated for error estimation
        calc = Calculator(w_data, n_bootstrap=50, step_size=1) 
        
        print(f"\nCalculating Effective Mass for R={R_target}")
        print("-" * 40)
        print(f"{'T':<5} | {'m_eff(T)':<12} | {'Error':<12}")
        print("-" * 40)
        
        # 3. Calculate and Print results
        # We loop through T values to see the effective mass plateau
        # m_eff(T) = ln( W(R,T) / W(R, T+1) )
        for t in range(1, 20):
            try:
                # The get_variable method handles calculation, caching, and bootstrapping
                var = calc.get_variable("effective_mass", R=R_target, T=t)
                
                val = var.get()
                err = var.err()
                
                # Check for valid numerical results
                if val is not None and not np.isnan(val):
                    # Format: T | Value | Error
                    print(f"{t:<5} | {val:.6f}     | {err:.6f}")
                else:
                    # If we hit NaNs (often due to negative loops at large T), stop output
                    pass
            except Exception:
                # Skip T values where calculation fails (e.g. missing data)
                continue
                
    except FileNotFoundError:
        print(f"Error: The directory '{file_path}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
