#!/usr/bin/env python3
"""
Read and store parameters from input.yaml into Python dataclasses.
"""

from dataclasses import dataclass
from typing import List, Any
import yaml
import argparse

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
class GradientFlowParams:
    enabled: bool
    integrator: str
    dt: float
    t_values: List[float]
    measure_energy_clover: bool
    measure_wilson_loop_temporal: bool
    measure_wilson_loop_mu_nu: bool
    extract_t0: bool
    t0_target: float
    obs_filename: str
    W_temp_filename: str
    W_mu_nu_filename: str
    t0_filename: str


@dataclass
class GaugeObservableParams:
    measurement_interval: int
    measure_plaquette: bool
    measure_wilson_loop_temporal: bool
    measure_wilson_loop_mu_nu: bool
    measure_retrace_U: bool
    W_temp_L_T_pairs: List[Any]
    W_mu_nu_pairs: List[Any]
    W_Lmu_Lnu_pairs: List[Any]
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

    metro_raw = dict(data["MetropolisParams"])
    gauge_raw = data["GaugeObservableParams"]

    metro_raw.setdefault("Ndims", 4)
    metro_raw.setdefault("Nd", 4)
    metro_raw.setdefault("Nc", 2)
    metro = MetropolisParams(**metro_raw)

    # Ensure the list entries are tuples (dataclass type annotation)
    # and expand string ranges (e.g. "1:4") into list of ints and take Cartesian product
    def _parse_range(val: Any) -> list[int]:
        if isinstance(val, int):
            return [val]
        elif isinstance(val, str) and ":" in val:
            parts = val.split(":")
            if len(parts) == 2:
                start, end = int(parts[0]), int(parts[1])
                return list(range(start, end + 1))
        # fallback, just try to parse as int
        return [int(val)]

    def _pairs(key: str) -> list[tuple[int, int]]:
        raw_list = gauge_raw.get(key, [])
        expanded = []
        for pair in raw_list:
            if len(pair) != 2:
                continue
            l_vals = _parse_range(pair[0])
            t_vals = _parse_range(pair[1])
            for l in l_vals:
                for t in t_vals:
                    expanded.append((l, t))
        return expanded

    gauge = GaugeObservableParams(
        measurement_interval=gauge_raw.get("measurement_interval", 1),
        measure_plaquette=gauge_raw.get("measure_plaquette", True),
        measure_wilson_loop_temporal=gauge_raw.get("measure_wilson_loop_temporal", False),
        measure_wilson_loop_mu_nu=gauge_raw.get("measure_wilson_loop_mu_nu", False),
        measure_retrace_U=gauge_raw.get("measure_retrace_U", False),
        W_temp_L_T_pairs=_pairs("W_temp_L_T_pairs"),
        W_mu_nu_pairs=_pairs("W_mu_nu_pairs"),
        W_Lmu_Lnu_pairs=_pairs("W_Lmu_Lnu_pairs"),
        plaquette_filename=gauge_raw.get("plaquette_filename", ""),
        W_temp_filename=gauge_raw.get("W_temp_filename", ""),
        W_mu_nu_filename=gauge_raw.get("W_mu_nu_filename", ""),
        RetraceU_filename=gauge_raw.get("RetraceU_filename", ""),
        write_to_file=gauge_raw.get("write_to_file", True),
    )

    return metro, gauge


def load_gradient_flow_params(yaml_path: str) -> GradientFlowParams:
    """Read optional GradientFlowParams from YAML with disabled defaults."""
    with open(yaml_path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)

    raw = dict(data.get("GradientFlowParams") or {})
    return GradientFlowParams(
        enabled=raw.get("enabled", False),
        integrator=raw.get("integrator", ""),
        dt=float(raw.get("dt", 0.0)),
        t_values=[float(value) for value in raw.get("t_values", [])],
        measure_energy_clover=raw.get("measure_energy_clover", False),
        measure_wilson_loop_temporal=raw.get("measure_wilson_loop_temporal", False),
        measure_wilson_loop_mu_nu=raw.get("measure_wilson_loop_mu_nu", False),
        extract_t0=raw.get("extract_t0", False),
        t0_target=float(raw.get("t0_target", 0.3)),
        obs_filename=raw.get("obs_filename", "gradient_flow_obs.dat"),
        W_temp_filename=raw.get("W_temp_filename", "gradient_flow_wtemp.dat"),
        W_mu_nu_filename=raw.get("W_mu_nu_filename", "gradient_flow_w_mu_nu.dat"),
        t0_filename=raw.get("t0_filename", "gradient_flow_t0.dat"),
    )


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
