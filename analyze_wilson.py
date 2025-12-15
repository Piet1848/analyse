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
    
    # New field for autocorrelation
    tau_int: float | None = None
    
    dataset: str = "Unknown"
    beta: float | None = None
    epsilon1: float | None = None
    epsilon2: float | None = None
    lattice_extents: Tuple[int, int, int, int] | None = None


def calculate_tau_int(series: np.ndarray, S: float = 1.5) -> float:
    """
    Calculate the integrated autocorrelation time tau_int using 
    automatic windowing (W >= S * tau).
    """
    n = len(series)
    if n < 100: return 0.5 # Not enough data
    
    # Subtract mean
    series = series - np.mean(series)
    var = np.var(series)
    if var == 0: return 0.5

    # Calculate autocorrelation function rho(t) using FFT
    ft = np.fft.rfft(series)
    spec = np.abs(ft)**2
    acf = np.fft.irfft(spec)
    # Only keep first half
    acf = acf[:n//2] 
    
    # Normalize
    # We use simple N normalization for the FFT estimate
    rho = (acf / np.arange(n, n - len(acf), -1)) / var 

    # Integrate rho to find tau
    tau = 0.5
    for t in range(1, len(rho)):
        tau += rho[t]
        # Automatic Windowing condition
        if t >= S * tau:
            return tau
            
    return tau


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
        if T >= t_min and val > 0: 
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
            
            popt, pcov = curve_fit(
                exponential_ansatz, 
                ts, 
                ws, 
                p0=[p0_C, p0_V],
                bounds=([-np.inf, -10.0], [np.inf, np.inf]) 
            )
            
            V_fit = popt[1]
            V_err = np.sqrt(np.diag(pcov))[1]
            
            potentials[L] = V_fit
            errors[L] = V_err
        except (RuntimeError, ValueError, Warning):
            pass
            
    return potentials, errors


def fit_sommer_parameter(
    potentials: Dict[int, float], 
    errors: Dict[int, float] = None,
    target_force_r2: float = 1.65
) -> Tuple[float | None, Dict[str, float] | None]:
    
    Ls = sorted(potentials.keys())
    if len(Ls) < 3:
        return None, None
        
    rs = np.array(Ls, dtype=float)
    Vs = np.array([potentials[L] for L in Ls])
    
    sigma = None
    if errors is not None:
        sigma = np.array([errors.get(L, 1.0) for L in Ls])
        avg_err = np.mean([e for e in sigma if e > 0]) if np.any(sigma > 0) else 1.0
        sigma = np.where(sigma <= 0, avg_err * 10, sigma)

    try:
        popt, pcov = curve_fit(
            cornell_potential_ansatz, 
            rs, 
            Vs, 
            p0=[0.0, 0.1, 0.26],
            sigma=sigma,             
            absolute_sigma=(sigma is not None)
        )
        
        A, sigma_val, B = popt
        
        numerator = target_force_r2 - B
        if numerator < 0 or sigma_val <= 0:
            return None, {'A': A, 'sigma': sigma_val, 'B': B}
            
        r0 = np.sqrt(numerator / sigma_val)
        
        return r0, {'A': A, 'sigma': sigma_val, 'B': B}
        
    except (RuntimeError, ValueError):
        return None, None