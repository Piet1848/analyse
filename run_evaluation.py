import data_organizer as do
import analyze_wilson as wilson
import argparse
import math
import numpy as np

def sommer_parameter(data: do.ExperimentData, sommer_target: float = 1.65):
    # Depending on how your data_organizer works, you might need to adjust this access
    # Assuming .get("W_temp") returns a list of ObservableData sets
    fileData = data.data.get("W_temp")[0] 
    
    # 1. Create Bootstrap samples
    fileData_bootstrap = fileData.get_bootstrap(n_bootstrap=200, seed=42)

    sommer_parameters = []
    lattice_spacings = []

    print(f"Starting bootstrap analysis on {len(fileData_bootstrap)} samples...")

    for i, fd in enumerate(fileData_bootstrap):
        # 2. Load records from the bootstrap sample
        records = wilson.load_w_temp(fd.observables)
        
        # 3. Average the Wilson loops W(R,T)
        averages = wilson.average_wilson_loops(records)

        # 4. Fit V(R) from the time dependence (NEW FUNCTION)
        # Returns potentials dict and errors dict
        potentials, _ = wilson.fit_potential_from_time(averages, t_min=2)

        # 5. Fit the Sommer parameter r0/a (NEW FUNCTION)
        # Returns r0 and the dictionary of Cornell parameters
        r0_over_a, _ = wilson.fit_sommer_parameter(potentials, target_force_r2=sommer_target)

        if r0_over_a and r0_over_a > 0:
            lattice_spacing = 0.5 / r0_over_a # Assuming r0_phys = 0.5 fm
            sommer_parameters.append(r0_over_a)
            lattice_spacings.append(lattice_spacing)
    
    if not lattice_spacings:
        print("Could not determine Sommer parameter for any bootstrap sample.")
        return None, None

    # Compute mean and stddev for Lattice Spacing
    mean_a = np.mean(lattice_spacings)
    stddev_a = np.std(lattice_spacings, ddof=1) # ddof=1 for sample stddev

    # Compute mean and stddev for r0/a
    mean_r0 = np.mean(sommer_parameters)
    stddev_r0 = np.std(sommer_parameters, ddof=1)

    print(f"Results over {len(lattice_spacings)} successful bootstrap samples:")
    print(f"  r0/a:               {mean_r0:.5f} ± {stddev_r0:.5f}")
    print(f"  Lattice spacing a:  {mean_a:.5f} ± {stddev_a:.5f} fm")

    return mean_a, stddev_a


if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument("path", type=str, help="Path to dataset folder")
    parsed_args = args.parse_args()
    
    path = parsed_args.path
    # path = "../data/20251205/3" 
    
    try:
        exp_data = do.ExperimentData(path)
        sommer_parameter(exp_data)
    except Exception as e:
        print(f"An error occurred: {e}")