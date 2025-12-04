#!/usr/bin/env python3
"""
Analyze temporal Wilson loops to extract V(R), the force F(r), and the Sommer
parameter r0/a for one or multiple dataset folders.
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

from load_input_yaml import load_params
import data_organizer as do


@dataclass
class WilsonLoopRecord:
    step: int
    L: int
    T: int
    value: float


@dataclass
class DatasetResult:
    averages: Dict[Tuple[int, int], float]
    potentials: Dict[int, Dict[int, float]]
    effective_potential: Dict[int, float]
    effective_sources: Dict[int, int]
    forces: List[Tuple[float, float]]
    r0_over_a: float | None
    lattice_spacing_fm: float | None
    physical_extent_fm: Tuple[float, float, float, float] | None
    physical_volume_fm4: float | None


def load_w_temp(WLoop: list[do.ObservableData]) -> List[WilsonLoopRecord]:
    """Convert ObservableData to a list of WilsonLoopRecord entries."""
    records: List[WilsonLoopRecord] = []

    steps  = [wl for wl in WLoop if wl.name == "# step"]
    Ls     = [wl for wl in WLoop if wl.name == "L"]
    Ts     = [wl for wl in WLoop if wl.name == "T"]
    values = [wl for wl in WLoop if wl.name == "W_temp"]

    for step, L, T, val in zip(steps[0].values, Ls[0].values, Ts[0].values, values[0].values):
        record = WilsonLoopRecord(
            step=int(step),
            L=int(L),
            T=int(T),
            value=float(val)
        )
        records.append(record)
    
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
