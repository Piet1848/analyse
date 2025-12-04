import data_organizer as do
import analyze_wilson as wilson
import argparse
import math

def sommer_parameter(data: do.ExperimentData, sommer_target: float = 1.65):
    fileData = data.data.get("W_temp")[0]
    fileData_bootstrap = fileData.get_bootstrap(n_bootstrap=100, seed=42)

    sommer_parameters = []
    for fd in fileData_bootstrap:
        records = wilson.load_w_temp(fd.observables)
        averages = wilson.average_wilson_loops(records)

        potentials = wilson.compute_potentials(averages)
        effective, sources = wilson.select_effective_potential(potentials)
        forces = wilson.compute_force(effective)
        r0_over_a = wilson.interpolate_sommer_parameter(forces, target=sommer_target)

        if r0_over_a and r0_over_a > 0:
            lattice_spacing = 0.5 / r0_over_a

        sommer_parameters.append(lattice_spacing)
    
    # Compute mean and stddev
    mean_a = sum(sommer_parameters) / len(sommer_parameters)
    stddev_a = math.sqrt(sum((x - mean_a) ** 2 for x in sommer_parameters) / (len(sommer_parameters) - 1))
    print(f"Sommer parameter results over {len(sommer_parameters)} bootstrap samples:")
    print(f"  Lattice spacing a (fm): {mean_a:.5f} ± {stddev_a:.5f}")

    return mean_a, stddev_a


if __name__ == "__main__":
    '''
    args = argparse.ArgumentParser()
    args.add_argument("path", type=str, help="Path to dataset folder")
    parsed_args = args.parse_args()
    path = parsed_args.path
    '''
    path = "../data/20251201/6"
    exp_data = do.ExperimentData(path)
    sommer_parameter(exp_data)