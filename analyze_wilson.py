#!/usr/bin/env python3
"""
Analyze temporal Wilson loops to extract V(R), the force F(r), and the Sommer
parameter r0/a using curve fitting.
"""

from __future__ import annotations

import math
import numpy as np
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple, Optional
from scipy.optimize import curve_fit

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
    effective_potential: Dict[int, float] 
    potential_errors: Dict[int, float]
    
    cornell_params: Dict[str, float] | None 
    
    r0_over_a: float | None
    lattice_spacing_fm: float | None
    physical_extent_fm: Tuple[float, float, float, float] | None
    physical_volume_fm4: float | None
    
    dataset: str = "Unknown"
    beta: float | None = None
    epsilon1: float | None = None
    epsilon2: float | None = None
    lattice_extents: Tuple[int, int, int, int] | None = None


def load_w_temp(WLoop: list[do.ObservableData]) -> List[WilsonLoopRecord]:
    """Convert ObservableData to a list of WilsonLoopRecord entries."""
    records: List[WilsonLoopRecord] = []
    steps  = [wl for wl in WLoop if wl.name == "# step"]
    Ls     = [wl for wl in WLoop if wl.name == "L"]
    Ts     = [wl for wl in WLoop if wl.name == "T"]
    values = [wl for wl in WLoop if wl.name == "W_temp"]

    if not steps or not Ls or not Ts or not values:
        return []

    for step, L, T, val in zip(steps[0].values, Ls[0].values, Ts[0].values, values[0].values):
        record = WilsonLoopRecord(step=int(step), L=int(L), T=int(T), value=float(val))
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


def exponential_ansatz(t, C, V):
    return C * np.exp(-V * t)

def cornell_potential_ansatz(r, A, sigma, B):
    return A + sigma * r - B / r

def fit_potential_from_time(
    avg: Dict[Tuple[int, int], float], 
    t_min: int = 2
) -> Tuple[Dict[int, float], Dict[int, float]]:
    # Group data by L
    data_by_L = defaultdict(list)
    for (L, T), val in avg.items():
        if T >= t_min and val > 0: # Filter early times and bad data
            data_by_L[L].append((T, val))
    
    potentials = {}
    errors = {}

    for L, points in data_by_L.items():
        if len(points) < 3:
            continue
            
        points.sort(key=lambda x: x[0])
        ts = np.array([p[0] for p in points])
        ws = np.array([p[1] for p in points])
        
        try:
            # Initial guess
            p0_V = -np.log(ws[1]/ws[0]) if ws[0] > 0 and ws[1] > 0 else 0.5
            p0_C = ws[0] * np.exp(p0_V * ts[0])
            
            popt, pcov = curve_fit(exponential_ansatz, ts, ws, p0=[p0_C, p0_V])
            
            V_fit = popt[1]
            V_err = np.sqrt(np.diag(pcov))[1]
            
            potentials[L] = V_fit
            errors[L] = V_err
        except (RuntimeError, ValueError, Warning):
            # Fallback or skip if fit fails
            pass
            
    return potentials, errors


def fit_sommer_parameter(
    potentials: Dict[int, float], 
    target_force_r2: float = 1.65
) -> Tuple[float | None, Dict[str, float] | None]:
    Ls = sorted(potentials.keys())
    if len(Ls) < 3:
        return None, None
        
    rs = np.array(Ls, dtype=float)
    Vs = np.array([potentials[L] for L in Ls])
    
    try:
        popt, _ = curve_fit(cornell_potential_ansatz, rs, Vs, p0=[0.0, 0.1, 0.26])
        A, sigma, B = popt
        
        numerator = target_force_r2 - B
        if numerator < 0 or sigma <= 0:
            return None, {'A': A, 'sigma': sigma, 'B': B}
            
        r0 = np.sqrt(numerator / sigma)
        
        return r0, {'A': A, 'sigma': sigma, 'B': B}
        
    except (RuntimeError, ValueError):
        return None, None

# --- Formatting ---
def format_dataset_result(result: DatasetResult) -> str:
    lines: List[str] = []
    header = f"Dataset: {result.dataset}"
    if result.beta is not None:
        header += f" (beta={result.beta})"
    lines.append(header)

    lines.append("  Fitted Potentials V(R) [from W(T) = C*exp(-VT)]:")
    for L in sorted(result.effective_potential):
        V = result.effective_potential[L]
        err = result.potential_errors.get(L, 0.0)
        lines.append(f"    R={L:>2d}: V={V:.6f} +/- {err:.6f}")

    if result.cornell_params:
        p = result.cornell_params
        lines.append("  Cornell Fit V(r) = A + sigma*r - B/r:")
        lines.append(f"    sigma={p['sigma']:.4f}, B={p['B']:.4f}, A={p['A']:.4f}")

    if result.r0_over_a is not None:
        lines.append(f"  Sommer parameter r0/a: {result.r0_over_a:.4f}")
        lines.append(f"  Lattice spacing a [fm] (assuming r0=0.5 fm): {result.lattice_spacing_fm:.4f}")
    else:
        lines.append("  Sommer parameter r0/a: Could not be determined.")
        
    return "\n".join(lines)