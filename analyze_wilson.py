#!/usr/bin/env python3
"""
Analyze temporal Wilson loops to extract V(R), the force F(r), and the Sommer
parameter r0/a for one or multiple dataset folders.
"""

from __future__ import annotations

import argparse
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

from load_input_yaml import load_params


@dataclass
class WilsonLoopRecord:
    step: int
    L: int
    T: int
    value: float


@dataclass
class DatasetResult:
    dataset: Path
    beta: float | None
    epsilon1: float | None
    epsilon2: float | None
    lattice_extents: Tuple[int, int, int, int] | None
    averages: Dict[Tuple[int, int], float]
    potentials: Dict[int, Dict[int, float]]
    effective_potential: Dict[int, float]
    effective_sources: Dict[int, int]
    forces: List[Tuple[float, float]]
    r0_over_a: float | None
    lattice_spacing_fm: float | None
    physical_extent_fm: Tuple[float, float, float, float] | None
    physical_volume_fm4: float | None


def read_w_temp(path: Path, start_step: int = 0) -> List[WilsonLoopRecord]:
    """Read W_temp.out as a list of WilsonLoopRecord entries."""
    records: List[WilsonLoopRecord] = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = [p.strip() for p in stripped.split(",")]
            if len(parts) != 4:
                raise ValueError(f"Unexpected line in {path}: {line}")
            step, L, T = map(int, parts[:3])
            if step < start_step:
                continue
            value = float(parts[3])
            records.append(WilsonLoopRecord(step=step, L=L, T=T, value=value))
    return records


def average_wilson_loops(records: Iterable[WilsonLoopRecord]) -> Dict[Tuple[int, int], float]:
    """Average W(L, T) over all measurements for each (L, T)."""
    sums: Dict[Tuple[int, int], float] = defaultdict(float)
    counts: Dict[Tuple[int, int], int] = defaultdict(int)
    for rec in records:
        key = (rec.L, rec.T)
        sums[key] += rec.value
        counts[key] += 1
    return {key: sums[key] / counts[key] for key in sums}


def compute_potentials(avg: Dict[Tuple[int, int], float]) -> Dict[int, Dict[int, float]]:
    """
    Compute V(R, T) = -ln(W(R, T+1) / W(R, T)) using averaged Wilson loops.
    Returns a mapping L -> {T: V(L, T)} for all available T pairs.
    """
    potentials: Dict[int, Dict[int, float]] = defaultdict(dict)
    for (L, T), w_rt in avg.items():
        next_key = (L, T + 1)
        if next_key not in avg:
            continue
        w_next = avg[next_key]
        if w_rt <= 0 or w_next <= 0:
            continue
        potentials[L][T] = -math.log(w_next / w_rt)
    return potentials


def select_effective_potential(
    potentials: Dict[int, Dict[int, float]]
) -> Tuple[Dict[int, float], Dict[int, int]]:
    """
    Pick an effective potential V(L) for each L using the largest available T
    (i.e., deepest in Euclidean time) as a simple plateau estimate.
    Returns the chosen V(L) and the corresponding T for transparency.
    """
    effective: Dict[int, float] = {}
    sources: Dict[int, int] = {}
    for L, t_map in potentials.items():
        if not t_map:
            continue
        t_star = max(t_map.keys())
        effective[L] = t_map[t_star]
        sources[L] = t_star
    return effective, sources


def compute_force(effective: Dict[int, float]) -> List[Tuple[float, float]]:
    """
    Compute a discrete force F(r) = dV/dr via nearest-neighbor differences.
    Returns a list of (r_midpoint, F) sorted by r_midpoint.
    """
    forces: List[Tuple[float, float]] = []
    if len(effective) < 2:
        return forces
    sorted_L = sorted(effective.keys())
    for left, right in zip(sorted_L[:-1], sorted_L[1:]):
        delta_r = right - left
        if delta_r == 0:
            continue
        r_mid = (left + right) / 2.0
        F = (effective[right] - effective[left]) / delta_r
        forces.append((r_mid, F))
    return forces


def interpolate_sommer_parameter(
    force_points: List[Tuple[float, float]], target: float = 1.65
) -> float | None:
    """
    Solve r^2 F(r) = target by linear interpolation between neighboring force
    points. Returns r0/a in lattice units, or None if no crossing is found.
    """
    if len(force_points) < 2:
        return None
    sorted_points = sorted(force_points, key=lambda x: x[0])
    for (r1, f1), (r2, f2) in zip(sorted_points[:-1], sorted_points[1:]):
        y1 = (r1**2) * f1
        y2 = (r2**2) * f2
        if (y1 - target) == 0:
            return r1
        if (y1 - target) * (y2 - target) > 0:
            continue
        if y2 == y1:
            continue
        t = (target - y1) / (y2 - y1)
        r0 = r1 + t * (r2 - r1)
        return r0
    return None


def analyze_dataset(dataset: Path, start_step: int = 0, sommer_target: float = 1.65) -> DatasetResult:
    """Full analysis pipeline for one dataset folder."""
    w_temp_path = dataset / "W_temp.out"
    if not w_temp_path.exists():
        raise FileNotFoundError(f"No W_temp.out found in {dataset}")

    records = read_w_temp(w_temp_path, start_step=start_step)
    averages = average_wilson_loops(records)
    potentials = compute_potentials(averages)
    effective, sources = select_effective_potential(potentials)
    forces = compute_force(effective)
    r0_over_a = interpolate_sommer_parameter(forces, target=sommer_target)
    lattice_spacing = None
    if r0_over_a and r0_over_a > 0:
        lattice_spacing = 0.5 / r0_over_a

    beta = None
    epsilon1 = None
    epsilon2 = None
    lattice_extents = None
    physical_extent_fm = None
    physical_volume_fm4 = None
    input_path = dataset / "input.yaml"
    if input_path.exists():
        metro, _ = load_params(str(input_path))
        beta = metro.beta
        epsilon1 = metro.epsilon1
        epsilon2 = metro.epsilon2
        lattice_extents = (metro.L0, metro.L1, metro.L2, metro.L3)
        if lattice_spacing:
            Ls = [d * lattice_spacing for d in lattice_extents]
            physical_extent_fm = tuple(Ls)  # type: ignore[assignment]
            physical_volume_fm4 = Ls[0] * Ls[1] * Ls[2] * Ls[3]

    return DatasetResult(
        dataset=dataset,
        beta=beta,
        epsilon1=epsilon1,
        epsilon2=epsilon2,
        lattice_extents=lattice_extents,
        averages=averages,
        potentials=potentials,
        effective_potential=effective,
        effective_sources=sources,
        forces=forces,
        r0_over_a=r0_over_a,
        lattice_spacing_fm=lattice_spacing,
        physical_extent_fm=physical_extent_fm,
        physical_volume_fm4=physical_volume_fm4,
    )


def discover_datasets(base_path: Path) -> List[Path]:
    """Return all subdirectories under base_path that contain W_temp.out."""
    if (base_path / "W_temp.out").exists():
        return [base_path]
    datasets = [
        p for p in base_path.iterdir() if p.is_dir() and (p / "W_temp.out").exists()
    ]
    return sorted(datasets)


def format_dataset_result(result: DatasetResult) -> str:
    lines: List[str] = []
    header = f"Dataset: {result.dataset}"
    if result.beta is not None:
        header += f" (beta={result.beta})"
    if result.epsilon1 is not None and result.epsilon2 is not None:
        header += f" eps1={result.epsilon1}, eps2={result.epsilon2}"
    lines.append(header)

    if result.lattice_extents:
        L0, L1, L2, L3 = result.lattice_extents
        lines.append(f"  Lattice sizes: L0={L0}, L1={L1}, L2={L2}, L3={L3}")

    lines.append("  Effective potentials V(R) [using largest available T]:")
    for L in sorted(result.effective_potential):
        V = result.effective_potential[L]
        t_used = result.effective_sources.get(L, "?")
        lines.append(f"    R={L:>2d}: V={V:.6f}  (from T={t_used}->T+1)")

    if result.forces:
        lines.append("  Forces F(r)=dV/dr and r^2 F(r):")
        for r_mid, F in result.forces:
            lines.append(f"    r={r_mid:4.1f}: F={F:.6f}, r^2F={r_mid*r_mid*F:.6f}")
    else:
        lines.append("  Forces: not enough data to compute.")

    if result.r0_over_a is not None:
        lines.append(f"  Sommer parameter r0/a: {result.r0_over_a:.4f}")
        lines.append(f"  Lattice spacing a [fm] (assuming r0=0.5 fm): {result.lattice_spacing_fm:.4f}")
        if result.physical_extent_fm:
            Ls = result.physical_extent_fm
            lines.append(
                "  Physical extents L*a [fm]: "
                f"L0={Ls[0]:.3f}, L1={Ls[1]:.3f}, L2={Ls[2]:.3f}, L3={Ls[3]:.3f}"
            )
        if result.physical_volume_fm4:
            lines.append(f"  Physical volume (L*a)^4 [fm^4]: {result.physical_volume_fm4:.3f}")
    else:
        lines.append("  Sommer parameter r0/a: not found (no crossing of r^2 F = 1.65).")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Analyze temporal Wilson loops and extract V(R), F(r), and r0/a."
    )
    parser.add_argument(
        "path",
        type=Path,
        help="Dataset folder or base directory containing dataset subfolders.",
    )
    parser.add_argument(
        "--start-step",
        type=int,
        default=0,
        help="Ignore measurements before this Monte Carlo step (thermalization cut).",
    )
    args = parser.parse_args()

    datasets = discover_datasets(args.path)
    if not datasets:
        raise SystemExit(f"No datasets with W_temp.out found under {args.path}")

    for dataset in datasets:
        result = analyze_dataset(dataset, start_step=args.start_step)
        print(format_dataset_result(result))
        print()


if __name__ == "__main__":
    main()
