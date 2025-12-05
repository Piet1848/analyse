#!/usr/bin/env python3
"""
Read and store parameters from input.yaml into Python dataclasses.
"""

from dataclasses import dataclass
from typing import List, Tuple
import yaml
import argparse


# ---------- Dataclass definitions ----------

@dataclass
class MetropolisParams:
    Ndims: int
    Nd: int
    Nc: int
    L0: int
    L1: int
    L2: int
    L3: int
    nHits: int
    nSweep: int
    seed: int
    beta: float
    delta: float
    epsilon1: float
    epsilon2: float


@dataclass
class GaugeObservableParams:
    measurement_interval: int
    measure_plaquette: bool
    measure_wilson_loop_temporal: bool
    measure_wilson_loop_mu_nu: bool
    measure_retrace_U: bool
    W_temp_L_T_pairs: List[Tuple[int, int]]
    W_mu_nu_pairs: List[Tuple[int, int]]
    W_Lmu_Lnu_pairs: List[Tuple[int, int]]
    plaquette_filename: str
    W_temp_filename: str
    W_mu_nu_filename: str
    RetraceU_filename: str
    write_to_file: bool

    # str
    def __str__(self) -> str:
        lines = [
            f"measurement_interval: {self.measurement_interval}",
            f"measure_plaquette: {self.measure_plaquette}",
            f"measure_wilson_loop_temporal: {self.measure_wilson_loop_temporal}",
            f"measure_wilson_loop_mu_nu: {self.measure_wilson_loop_mu_nu}",
            f"measure_retrace_U: {self.measure_retrace_U}",
            f"W_temp_L_T_pairs: {self.W_temp_L_T_pairs}",
            f"W_mu_nu_pairs: {self.W_mu_nu_pairs}",
            f"W_Lmu_Lnu_pairs: {self.W_Lmu_Lnu_pairs}",
            f"plaquette_filename: {self.plaquette_filename}",
            f"W_temp_filename: {self.W_temp_filename}",
            f"W_mu_nu_filename: {self.W_mu_nu_filename}",
            f"RetraceU_filename: {self.RetraceU_filename}",
            f"write_to_file: {self.write_to_file}",
        ]
        return "\n".join(lines)


# ---------- Load function ----------

def load_params(yaml_path: str) -> tuple[MetropolisParams, GaugeObservableParams]:
    """Read YAML file and return both parameter blocks as dataclasses."""
    with open(yaml_path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)

    metro_raw = data["MetropolisParams"]
    gauge_raw = data["GaugeObservableParams"]

    metro = MetropolisParams(**metro_raw)

    # Ensure the list entries are tuples (dataclass type annotation)
    def _pairs(key: str) -> list[tuple[int, int]]:
        return [tuple(pair) for pair in gauge_raw.get(key, [])]

    gauge = GaugeObservableParams(
        measurement_interval=gauge_raw["measurement_interval"],
        measure_plaquette=gauge_raw["measure_plaquette"],
        measure_wilson_loop_temporal=gauge_raw["measure_wilson_loop_temporal"],
        measure_wilson_loop_mu_nu=gauge_raw["measure_wilson_loop_mu_nu"],
        measure_retrace_U=gauge_raw["measure_retrace_U"],
        W_temp_L_T_pairs=_pairs("W_temp_L_T_pairs"),
        W_mu_nu_pairs=_pairs("W_mu_nu_pairs"),
        W_Lmu_Lnu_pairs=_pairs("W_Lmu_Lnu_pairs"),
        plaquette_filename=gauge_raw["plaquette_filename"],
        W_temp_filename=gauge_raw["W_temp_filename"],
        W_mu_nu_filename=gauge_raw["W_mu_nu_filename"],
        RetraceU_filename=gauge_raw["RetraceU_filename"],
        write_to_file=gauge_raw["write_to_file"],
    )

    return metro, gauge


if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument("path", type=str, help="Path to dataset folder")
    parsed_args = args.parse_args()
    path = parsed_args.path
    # Example usage
    metro_params, gauge_params = load_params(path + "/input.yaml")
    print("Metropolis Parameters:")
    print(metro_params)
    print("\nGauge Observable Parameters:")
    print(gauge_params)