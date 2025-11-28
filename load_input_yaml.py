#!/usr/bin/env python3
"""
Read and store parameters from input.yaml into Python dataclasses.
"""

from dataclasses import dataclass
from typing import List, Tuple
import yaml


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

