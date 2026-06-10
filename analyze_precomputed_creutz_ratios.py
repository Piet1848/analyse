#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline

from calculator import creutz_chi_from_wilson
from finalized_analysis_helpers import _get_plotly, _write_figure_html, save_matplotlib_figure


STATS_FILENAME = "wilson_loop_stats.dat"


@dataclass(frozen=True)
class WilsonLoopStat:
    flow_time: float
    L: int
    T: int
    mean: float
    error: float | None
    bootstrap_path: Path | None
    n_measurements: int | None = None


@dataclass(frozen=True)
class FitScaleConfig:
    name: str
    x_key: str
    y_key: str
    err_key: str
    scale_factor: float
    x_denominator: float
    x_title: str
    y_title: str
    x_label: str
    y_label: str


def filename_token(value: float | str) -> str:
    return str(value).replace("-", "m").replace(".", "p")


def flow_time_label(flow_time: float) -> str:
    return f"t_lat={float(flow_time):g}"


def flow_time_title(flow_time: float) -> str:
    return f"t_lat = {float(flow_time):g}"


def json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def save_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, default=json_default)
        handle.write("\n")


def discover_analysis_dirs(paths: list[Path]) -> list[Path]:
    seen: set[Path] = set()
    result: list[Path] = []

    for raw_path in paths:
        path = raw_path.expanduser().resolve()
        candidates: list[Path] = []
        if path.is_file() and path.name == STATS_FILENAME:
            candidates = [path.parent]
        elif (path / STATS_FILENAME).exists():
            candidates = [path]
        elif path.is_dir():
            candidates = sorted(item.parent for item in path.rglob(STATS_FILENAME))
        else:
            raise FileNotFoundError(f"No such input path: {path}")

        for candidate in candidates:
            if candidate not in seen:
                seen.add(candidate)
                result.append(candidate)

    return result


def parse_float(value: str) -> float | None:
    if value.lower() in {"", "none", "nan"}:
        return None
    return float(value)


def load_wilson_loop_stats(analysis_dir: Path) -> list[WilsonLoopStat]:
    table_path = analysis_dir / STATS_FILENAME
    if not table_path.exists():
        raise FileNotFoundError(f"Missing {STATS_FILENAME}: {table_path}")

    header: list[str] | None = None
    rows: list[WilsonLoopStat] = []
    with table_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                tokens = stripped[1:].strip().split()
                if tokens and tokens[:3] == ["t_over_a2", "L", "T"]:
                    header = tokens
                continue

            parts = stripped.split()
            if header is None:
                raise ValueError(f"{table_path}:{line_number}: missing header line")
            if len(parts) != len(header):
                raise ValueError(
                    f"{table_path}:{line_number}: expected {len(header)} columns, got {len(parts)}"
                )
            row = dict(zip(header, parts, strict=True))
            raw_bootstrap_path = row.get("bootstrap_path")
            bootstrap_path = None
            if raw_bootstrap_path:
                candidate = Path(raw_bootstrap_path)
                bootstrap_path = candidate if candidate.is_absolute() else analysis_dir / candidate

            rows.append(
                WilsonLoopStat(
                    flow_time=float(row["t_over_a2"]),
                    L=int(row["L"]),
                    T=int(row["T"]),
                    mean=float(row["mean"]),
                    error=parse_float(row.get("bootstrap_error", "nan")),
                    bootstrap_path=bootstrap_path,
                    n_measurements=int(row["n_measurements"]) if row.get("n_measurements") else None,
                )
            )

    return rows


def positive_log_ratio_error(
    w_rt: WilsonLoopStat,
    w_rt1: WilsonLoopStat,
    w_r1t: WilsonLoopStat,
    w_r1t1: WilsonLoopStat,
) -> float | None:
    values = [w_rt.mean, w_rt1.mean, w_r1t.mean, w_r1t1.mean]
    errors = [w_rt.error, w_rt1.error, w_r1t.error, w_r1t1.error]
    if any(value <= 0 or not np.isfinite(value) for value in values):
        return None
    if any(error is None or not np.isfinite(error) for error in errors):
        return None
    variance = sum((float(error) / value) ** 2 for value, error in zip(values, errors, strict=True))
    return float(math.sqrt(variance))


def load_bootstrap(stat: WilsonLoopStat) -> np.ndarray | None:
    if stat.bootstrap_path is None or not stat.bootstrap_path.exists():
        return None
    arr = np.asarray(np.load(stat.bootstrap_path), dtype=float)
    if arr.ndim != 1 or arr.size == 0:
        return None
    return arr


def bootstrap_creutz(
    w_rt: WilsonLoopStat,
    w_rt1: WilsonLoopStat,
    w_r1t: WilsonLoopStat,
    w_r1t1: WilsonLoopStat,
) -> tuple[np.ndarray | None, int]:
    samples = [load_bootstrap(stat) for stat in (w_rt, w_rt1, w_r1t, w_r1t1)]
    if any(sample is None for sample in samples):
        return None, 0
    arrays = [np.asarray(sample, dtype=float) for sample in samples if sample is not None]
    n_samples = min(arr.size for arr in arrays)
    if n_samples == 0:
        return None, 0

    w_rt_boot, w_rt1_boot, w_r1t_boot, w_r1t1_boot = [arr[:n_samples] for arr in arrays]
    chi_boot = np.asarray(
        creutz_chi_from_wilson(w_rt_boot, w_rt1_boot, w_r1t_boot, w_r1t1_boot),
        dtype=float,
    )
    finite = int(np.count_nonzero(np.isfinite(chi_boot)))
    if finite == 0:
        return None, 0
    return chi_boot, finite


def relative_to(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def compute_creutz_records(
    analysis_dir: Path,
    output_dir: Path,
    rows: list[WilsonLoopStat],
    *,
    t0_over_a2: float | None,
    r0_over_a: float | None,
    use_bootstrap: bool,
    diagonal_only: bool,
) -> list[dict[str, Any]]:
    by_flow_pair = {(row.flow_time, row.L, row.T): row for row in rows}
    flow_times = sorted({row.flow_time for row in rows})
    records: list[dict[str, Any]] = []
    bootstrap_root = output_dir / "bootstrap"

    for flow_time in flow_times:
        pairs = {(row.L, row.T) for row in rows if row.flow_time == flow_time}
        for l_size, t_size in sorted(pairs):
            required = [
                (l_size, t_size),
                (l_size, t_size + 1),
                (l_size + 1, t_size),
                (l_size + 1, t_size + 1),
            ]
            if not all(pair in pairs for pair in required):
                continue

            r_mid = float(l_size) + 0.5
            t_mid = float(t_size) + 0.5
            if diagonal_only and not np.isclose(r_mid, t_mid):
                continue

            w_rt = by_flow_pair[(flow_time, l_size, t_size)]
            w_rt1 = by_flow_pair[(flow_time, l_size, t_size + 1)]
            w_r1t = by_flow_pair[(flow_time, l_size + 1, t_size)]
            w_r1t1 = by_flow_pair[(flow_time, l_size + 1, t_size + 1)]

            chi = creutz_chi_from_wilson(w_rt.mean, w_rt1.mean, w_r1t.mean, w_r1t1.mean)
            if chi is None or not np.isfinite(chi):
                continue

            chi_boot: np.ndarray | None = None
            n_bootstrap_valid = 0
            if use_bootstrap:
                chi_boot, n_bootstrap_valid = bootstrap_creutz(w_rt, w_rt1, w_r1t, w_r1t1)

            chi_err = None
            bootstrap_path = None
            if chi_boot is not None:
                chi_err = float(np.nanstd(chi_boot))
                flow_dir = bootstrap_root / f"t_over_a2_{filename_token(f'{flow_time:g}')}"
                flow_dir.mkdir(parents=True, exist_ok=True)
                bootstrap_path = flow_dir / (
                    f"chi_R_{filename_token(f'{r_mid:g}')}_T_{filename_token(f'{t_mid:g}')}.npy"
                )
                np.save(bootstrap_path, chi_boot)
            else:
                chi_err = positive_log_ratio_error(w_rt, w_rt1, w_r1t, w_r1t1)

            record: dict[str, Any] = {
                "t_over_a2": float(flow_time),
                "t_lat": float(flow_time),
                "eight_t_over_a2": float(8.0 * flow_time),
                "R": int(l_size),
                "T": int(t_size),
                "R_mid": float(r_mid),
                "T_mid": float(t_mid),
                "chi": float(chi),
                "chi_err": chi_err,
                "n_bootstrap_valid": int(n_bootstrap_valid),
                "chi_bootstrap_path": relative_to(bootstrap_path, output_dir) if bootstrap_path else None,
                "source_analysis_dir": str(analysis_dir),
            }

            if t0_over_a2 is not None:
                scale = float(8.0 * t0_over_a2)
                record["t0_over_a2"] = float(t0_over_a2)
                record["r_hat_t0"] = float(r_mid / math.sqrt(scale))
                record["t_hat_t0"] = float(t_mid / math.sqrt(scale))
                record["chi_hat_t0"] = float(chi * scale)
                record["chi_hat_t0_err"] = float(chi_err * scale) if chi_err is not None else None

            if r0_over_a is not None:
                scale = float(r0_over_a) ** 2
                record["r0_over_a"] = float(r0_over_a)
                record["r_hat_r0"] = float(r_mid / float(r0_over_a))
                record["t_hat_r0"] = float(t_mid / float(r0_over_a))
                record["chi_hat_r0"] = float(chi * scale)
                record["chi_hat_r0_err"] = float(chi_err * scale) if chi_err is not None else None

            records.append(record)

    return records


def write_table(path: Path, records: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "t_over_a2",
        "t_lat",
        "eight_t_over_a2",
        "R",
        "T",
        "R_mid",
        "T_mid",
        "chi",
        "chi_err",
        "r_hat_t0",
        "t_hat_t0",
        "chi_hat_t0",
        "chi_hat_t0_err",
        "r_hat_r0",
        "t_hat_r0",
        "chi_hat_r0",
        "chi_hat_r0_err",
        "n_bootstrap_valid",
        "chi_bootstrap_path",
    ]
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# " + " ".join(columns) + "\n")
        for record in records:
            values = []
            for column in columns:
                value = record.get(column)
                if value is None:
                    values.append("nan")
                elif isinstance(value, float):
                    values.append(f"{value:.12g}")
                else:
                    values.append(str(value))
            handle.write(" ".join(values) + "\n")


def diagonal_records(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        record
        for record in records
        if np.isclose(float(record["R_mid"]), float(record["T_mid"]))
    ]


def _finite_positive_sigma(sigma: np.ndarray) -> np.ndarray:
    sigma = np.asarray(sigma, dtype=float)
    positive = sigma[np.isfinite(sigma) & (sigma > 0)]
    fallback = float(np.median(positive)) if positive.size else 1.0
    result = sigma.copy()
    result[~np.isfinite(result) | (result <= 0)] = fallback
    return result


def _local_power_interpolation(
    x: np.ndarray,
    y: np.ndarray,
    x_grid: np.ndarray,
    powers: tuple[int, ...],
) -> tuple[np.ndarray, dict[str, Any]]:
    window_size = len(powers)
    if x.size < window_size:
        raise ValueError(f"need at least {window_size} points")

    order = np.argsort(x)
    x_sorted = np.asarray(x[order], dtype=float)
    y_sorted = np.asarray(y[order], dtype=float)
    if np.unique(x_sorted).size != x_sorted.size:
        raise ValueError("polynomial interpolation x values must be unique")

    window_coeffs: list[np.ndarray] = []
    window_ranges: list[tuple[float, float]] = []
    for start in range(0, len(x_sorted) - window_size + 1):
        x_window = x_sorted[start : start + window_size]
        y_window = y_sorted[start : start + window_size]
        design = np.column_stack([x_window**power for power in powers])
        coeffs = np.linalg.solve(design, y_window)
        window_coeffs.append(coeffs)
        window_ranges.append((float(x_window[0]), float(x_window[-1])))

    curve = np.full_like(np.asarray(x_grid, dtype=float), np.nan, dtype=float)
    for index, x_value in enumerate(np.asarray(x_grid, dtype=float)):
        values: list[float] = []
        for coeffs, (x_min, x_max) in zip(window_coeffs, window_ranges, strict=True):
            if x_min <= x_value <= x_max:
                basis = np.asarray([x_value**power for power in powers], dtype=float)
                values.append(float(basis @ coeffs))
        if not values:
            nearest_index = int(
                np.argmin(
                    [
                        min(abs(x_value - x_min), abs(x_value - x_max))
                        for x_min, x_max in window_ranges
                    ]
                )
            )
            coeffs = window_coeffs[nearest_index]
            basis = np.asarray([x_value**power for power in powers], dtype=float)
            values.append(float(basis @ coeffs))
        curve[index] = float(np.mean(values))

    return curve, {
        "kind": "local_adjacent_power_interpolation",
        "powers": list(powers),
        "window_size": int(window_size),
        "n_windows": len(window_coeffs),
        "window_coefficients": [[float(value) for value in coeffs] for coeffs in window_coeffs],
    }


def _natural_spline_interpolation(
    x: np.ndarray,
    y: np.ndarray,
    x_grid: np.ndarray,
    *,
    remove_inverse_square: bool = False,
) -> tuple[np.ndarray, dict[str, Any]]:
    if x.size < 3:
        raise ValueError("need at least 3 points")

    order = np.argsort(x)
    x_sorted = x[order]
    y_sorted = y[order]
    if np.unique(x_sorted).size != x_sorted.size:
        raise ValueError("spline x values must be unique")

    if remove_inverse_square:
        y_sorted = y_sorted * x_sorted**2

    spline = CubicSpline(x_sorted, y_sorted, bc_type="natural")
    curve = np.asarray(spline(x_grid), dtype=float)
    if remove_inverse_square:
        curve = curve / (x_grid**2)
    return curve, {
        "kind": "natural_cubic_spline",
        "remove_inverse_square": bool(remove_inverse_square),
    }


def fit_model_curves(
    x: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray,
    x_grid: np.ndarray,
) -> dict[str, Any]:
    specs = [
        (
            "cspl_pol_m2",
            "CSpl * Pol_-2",
            lambda target_x: _natural_spline_interpolation(x, y, target_x, remove_inverse_square=True),
        ),
        (
            "cspl",
            "CSpl",
            lambda target_x: _natural_spline_interpolation(x, y, target_x, remove_inverse_square=False),
        ),
        (
            "pol_0_m2_m4",
            "Pol_0,-2,-4",
            lambda target_x: _local_power_interpolation(x, y, target_x, (0, -2, -4)),
        ),
        (
            "pol_2_1_0",
            "Pol_2,1,0",
            lambda target_x: _local_power_interpolation(x, y, target_x, (2, 1, 0)),
        ),
    ]

    sigma = _finite_positive_sigma(sigma)
    models: list[dict[str, Any]] = []
    failed: list[dict[str, str]] = []
    for name, label, factory in specs:
        try:
            curve, metadata = factory(x_grid)
            y_at_points, _metadata_at_points = factory(x)
        except Exception as exc:
            failed.append({"name": name, "label": label, "reason": str(exc)})
            continue
        if not np.any(np.isfinite(curve)) or not np.any(np.isfinite(y_at_points)):
            failed.append({"name": name, "label": label, "reason": "non-finite curve"})
            continue
        residuals = (np.asarray(y_at_points, dtype=float) - y) / sigma
        finite_residuals = residuals[np.isfinite(residuals)]
        if finite_residuals.size == 0:
            failed.append({"name": name, "label": label, "reason": "non-finite residuals"})
            continue
        chi2 = float(np.sum(finite_residuals**2))
        chi2_dof = float(chi2 / max(1, finite_residuals.size - 1))
        max_abs_pull = float(np.max(np.abs(finite_residuals)))
        models.append(
            {
                "name": name,
                "label": label,
                "y": curve,
                "y_at_points": np.asarray(y_at_points, dtype=float),
                "chi2": chi2,
                "chi2_dof": chi2_dof,
                "max_abs_pull": max_abs_pull,
                "used_in_average": True,
                "metadata": metadata,
            }
        )

    if not models:
        raise ValueError("all fit models failed")

    matrix = np.asarray([np.asarray(model["y"], dtype=float) for model in models], dtype=float)
    average = np.nanmean(matrix, axis=0)
    systematic = (
        0.5 * (np.nanmax(matrix, axis=0) - np.nanmin(matrix, axis=0))
        if matrix.shape[0] > 1
        else np.zeros_like(average)
    )
    return {
        "models": models,
        "failed_models": failed,
        "average": average,
        "systematic": systematic,
    }


def filtered_fit_rows(
    rows: list[dict[str, Any]],
    *,
    x_key: str,
    y_key: str,
    err_key: str,
    max_rel_err: float | None,
) -> list[dict[str, Any]]:
    selected: list[dict[str, Any]] = []
    for row in rows:
        try:
            x_val = float(row[x_key])
            y_val = float(row[y_key])
            err = row.get(err_key)
            err_val = float(err) if err is not None else float("nan")
        except (KeyError, TypeError, ValueError):
            continue
        if not (np.isfinite(x_val) and np.isfinite(y_val) and np.isfinite(err_val)):
            continue
        if x_val <= 0 or y_val <= 0 or err_val <= 0:
            continue
        if max_rel_err is not None and abs(err_val / y_val) > max_rel_err:
            continue
        selected.append(row)
    return sorted(selected, key=lambda item: float(item[x_key]))


def load_dimensionless_bootstrap_matrix(
    rows: list[dict[str, Any]],
    output_dir: Path,
    *,
    scale_factor: float,
) -> np.ndarray | None:
    arrays: list[np.ndarray] = []
    for row in rows:
        raw_path = row.get("chi_bootstrap_path")
        if raw_path in {None, "", "nan"}:
            return None
        path = Path(str(raw_path))
        path = path if path.is_absolute() else output_dir / path
        if not path.exists():
            return None
        arr = np.asarray(np.load(path), dtype=float)
        if arr.ndim != 1 or arr.size == 0:
            return None
        arrays.append(arr * scale_factor)

    n_samples = min(arr.size for arr in arrays) if arrays else 0
    if n_samples <= 0:
        return None
    return np.vstack([arr[:n_samples] for arr in arrays])


def statistical_fit_band(
    fit_rows: list[dict[str, Any]],
    output_dir: Path,
    config: FitScaleConfig,
    x_grid: np.ndarray,
    *,
    max_samples: int,
    seed: int,
) -> tuple[np.ndarray, int, str]:
    x = np.asarray([float(row[config.x_key]) for row in fit_rows], dtype=float)
    y = np.asarray([float(row[config.y_key]) for row in fit_rows], dtype=float)
    sigma = np.asarray([float(row[config.err_key]) for row in fit_rows], dtype=float)

    bootstrap_matrix = load_dimensionless_bootstrap_matrix(
        fit_rows,
        output_dir,
        scale_factor=config.scale_factor,
    )
    source = "wilson_loop_bootstrap"
    if bootstrap_matrix is None:
        rng = np.random.default_rng(seed)
        n_generated = max(1, int(max_samples))
        bootstrap_matrix = rng.normal(y[:, None], sigma[:, None], size=(x.size, n_generated))
        source = "gaussian_from_errors"

    n_available = int(bootstrap_matrix.shape[1])
    n_samples = min(max(1, int(max_samples)), n_available)
    if n_samples < n_available:
        indices = np.linspace(0, n_available - 1, n_samples, dtype=int)
    else:
        indices = np.arange(n_available, dtype=int)

    average_curves: list[np.ndarray] = []
    for index in indices:
        y_replica = np.asarray(bootstrap_matrix[:, index], dtype=float)
        mask = (
            np.isfinite(x)
            & np.isfinite(y_replica)
            & np.isfinite(sigma)
            & (x > 0)
            & (y_replica > 0)
            & (sigma > 0)
        )
        if np.count_nonzero(mask) < 3:
            continue
        try:
            fitted = fit_model_curves(
                x[mask],
                y_replica[mask],
                sigma[mask],
                x_grid,
            )
        except Exception:
            continue
        average_curves.append(np.asarray(fitted["average"], dtype=float))

    if not average_curves:
        return np.zeros_like(x_grid, dtype=float), 0, source

    return np.nanstd(np.asarray(average_curves, dtype=float), axis=0), len(average_curves), source


def smoothing_radius_in_plot_units(flow_time: float, config: FitScaleConfig) -> float | None:
    eight_t_over_a2 = 8.0 * float(flow_time)
    if eight_t_over_a2 <= 0 or config.x_denominator <= 0:
        return None
    return float(math.sqrt(eight_t_over_a2) / config.x_denominator)


def build_dimensionless_fit_for_flow(
    rows: list[dict[str, Any]],
    output_dir: Path,
    config: FitScaleConfig,
    *,
    max_rel_err: float | None,
    bootstrap_samples: int,
    seed: int,
) -> dict[str, Any] | None:
    all_rows = sorted(
        [row for row in rows if config.x_key in row and config.y_key in row],
        key=lambda item: float(item[config.x_key]),
    )
    fit_rows = filtered_fit_rows(
        all_rows,
        x_key=config.x_key,
        y_key=config.y_key,
        err_key=config.err_key,
        max_rel_err=max_rel_err,
    )
    if len(fit_rows) < 3:
        return None

    x = np.asarray([float(row[config.x_key]) for row in fit_rows], dtype=float)
    y = np.asarray([float(row[config.y_key]) for row in fit_rows], dtype=float)
    sigma = np.asarray([float(row[config.err_key]) for row in fit_rows], dtype=float)
    x_min = float(np.min(x))
    x_max = float(np.max(x))
    x_grid = np.linspace(x_min, x_max, 300)
    central = fit_model_curves(
        x,
        y,
        sigma,
        x_grid,
    )
    stat, n_stat_samples, stat_source = statistical_fit_band(
        fit_rows,
        output_dir,
        config,
        x_grid,
        max_samples=bootstrap_samples,
        seed=seed,
    )
    systematic = np.asarray(central["systematic"], dtype=float)
    total = np.sqrt(stat**2 + systematic**2)
    average = np.asarray(central["average"], dtype=float)
    denominator = np.where(np.abs(average) > 0, np.abs(average), np.nan)

    flow_time = float(fit_rows[0]["t_over_a2"])
    return {
        "scale": config.name,
        "flow_time": flow_time,
        "t_lat": flow_time,
        "eight_t_over_a2": float(8.0 * flow_time),
        "n_points_available": len(all_rows),
        "n_points_fit": len(fit_rows),
        "fit_options": {
            "max_rel_err": float(max_rel_err) if max_rel_err is not None else None,
            "bootstrap_samples_requested": int(bootstrap_samples),
            "statistical_source": stat_source,
            "statistical_samples_used": int(n_stat_samples),
        },
        "smoothing_radius_r_hat": smoothing_radius_in_plot_units(flow_time, config),
        "points": [
            {
                "R_mid": float(row["R_mid"]),
                "T_mid": float(row["T_mid"]),
                "x": float(row[config.x_key]),
                "y": float(row[config.y_key]),
                "err": float(row[config.err_key]),
                "used_in_fit": row in fit_rows,
            }
            for row in all_rows
            if row.get(config.err_key) is not None
        ],
        "grid": {
            "x": x_grid,
            "average": average,
            "stat_err": stat,
            "sys_err": systematic,
            "total_err": total,
            "stat_rel": stat / denominator,
            "sys_rel": systematic / denominator,
            "total_rel": total / denominator,
            "models": {
                model["name"]: {
                    "label": model["label"],
                    "y": np.asarray(model["y"], dtype=float),
                    "used_in_average": bool(model["used_in_average"]),
                    "chi2_dof": float(model["chi2_dof"]),
                    "max_abs_pull": float(model["max_abs_pull"]),
                    "metadata": model["metadata"],
                }
                for model in central["models"]
            },
        },
        "failed_models": central["failed_models"],
    }


def save_fit_png(path: Path, fit: dict[str, Any], config: FitScaleConfig) -> dict[str, Path]:
    grid = fit["grid"]
    x_grid = np.asarray(grid["x"], dtype=float)
    fig, (ax_fit, ax_rel) = plt.subplots(1, 2, figsize=(9.2, 3.4), dpi=160)

    points = [point for point in fit["points"] if point.get("used_in_fit")] or fit["points"]
    ax_fit.errorbar(
        [float(point["x"]) for point in points],
        [float(point["y"]) for point in points],
        yerr=[float(point["err"]) for point in points],
        fmt="o",
        color="black",
        ecolor="black",
        markersize=2.4,
        capsize=1.4,
        linewidth=0.7,
        label="data",
    )
    ax_fit.plot(x_grid, np.asarray(grid["average"], dtype=float), color="0.35", linewidth=1.2, label="Average")
    model_styles = {
        "cspl_pol_m2": "#1f77b4",
        "cspl": "#ff7f0e",
        "pol_0_m2_m4": "#2ca02c",
        "pol_2_1_0": "#d62728",
    }
    for name, model in grid["models"].items():
        ax_fit.plot(
            x_grid,
            np.asarray(model["y"], dtype=float),
            linewidth=0.9,
            color=model_styles.get(name),
            label=str(model["label"]),
        )
    smoothing_radius = fit.get("smoothing_radius_r_hat")
    if smoothing_radius is not None and np.isfinite(smoothing_radius):
        ax_fit.axvline(float(smoothing_radius), color="0.65", linewidth=0.8)

    ax_rel.plot(x_grid, np.asarray(grid["total_rel"], dtype=float), color="#1f77b4", linewidth=1.0, label="Delta_tot")
    ax_rel.plot(x_grid, np.asarray(grid["sys_rel"], dtype=float), color="#2ca02c", linewidth=1.0, label="Delta_sys")
    ax_rel.plot(x_grid, np.asarray(grid["stat_rel"], dtype=float), color="#d62728", linewidth=1.0, label="Delta_stat")

    title = flow_time_title(float(fit.get("t_lat", fit["flow_time"])))
    ax_fit.text(
        0.50,
        0.95,
        title,
        ha="center",
        va="top",
        transform=ax_fit.transAxes,
        bbox={"facecolor": "white", "edgecolor": "0.4", "linewidth": 0.6},
        fontsize=8,
    )
    ax_fit.set_xlabel(config.x_label)
    ax_fit.set_ylabel(config.y_label)
    ax_rel.set_xlabel(config.x_label)
    ax_rel.set_ylabel(r"$\Delta\hat\chi/\hat\chi$")
    ax_fit.tick_params(direction="in", top=True, right=True)
    ax_rel.tick_params(direction="in", top=True, right=True)
    ax_fit.legend(frameon=True, fontsize=7)
    ax_rel.legend(frameon=True, fontsize=7)
    ax_rel.set_ylim(bottom=0)
    fig.tight_layout()
    written = save_matplotlib_figure(fig, path)
    plt.close(fig)
    if not written:
        raise RuntimeError(f"No fit plot exports could be written for {path}.")
    return written


def save_fit_html(path: Path, fit: dict[str, Any], config: FitScaleConfig) -> None:
    go = _get_plotly()
    from plotly.subplots import make_subplots

    grid = fit["grid"]
    x_grid = [float(value) for value in grid["x"]]
    fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=(
            flow_time_title(float(fit.get("t_lat", fit["flow_time"]))),
            "Relative fit uncertainty",
        ),
    )
    points = [point for point in fit["points"] if point.get("used_in_fit")] or fit["points"]
    fig.add_trace(
        go.Scatter(
            x=[float(point["x"]) for point in points],
            y=[float(point["y"]) for point in points],
            error_y={
                "type": "data",
                "array": [float(point["err"]) for point in points],
                "visible": True,
            },
            mode="markers",
            name="data",
            marker={"color": "black", "size": 5},
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=x_grid,
            y=[float(value) for value in grid["average"]],
            mode="lines",
            name="Average",
            line={"color": "gray", "width": 2},
        ),
        row=1,
        col=1,
    )
    model_colors = {
        "cspl_pol_m2": "#1f77b4",
        "cspl": "#ff7f0e",
        "pol_0_m2_m4": "#2ca02c",
        "pol_2_1_0": "#d62728",
    }
    for name, model in grid["models"].items():
        fig.add_trace(
            go.Scatter(
                x=x_grid,
                y=[float(value) for value in model["y"]],
                mode="lines",
                name=str(model["label"]),
                line={
                    "color": model_colors.get(name),
                    "width": 1.3,
                },
            ),
            row=1,
            col=1,
        )

    smoothing_radius = fit.get("smoothing_radius_r_hat")
    if smoothing_radius is not None and np.isfinite(smoothing_radius):
        y_values = [float(point["y"]) for point in points if np.isfinite(float(point["y"]))]
        y_values.extend(float(value) for value in grid["average"] if np.isfinite(float(value)))
        if y_values:
            fig.add_trace(
                go.Scatter(
                    x=[float(smoothing_radius), float(smoothing_radius)],
                    y=[min(y_values), max(y_values)],
                    mode="lines",
                    showlegend=False,
                    line={"color": "rgba(120,120,120,0.65)", "width": 1},
                    hoverinfo="skip",
                ),
                row=1,
                col=1,
            )

    for key, label, color in (
        ("total_rel", "Delta_tot", "#1f77b4"),
        ("sys_rel", "Delta_sys", "#2ca02c"),
        ("stat_rel", "Delta_stat", "#d62728"),
    ):
        fig.add_trace(
            go.Scatter(
                x=x_grid,
                y=[float(value) for value in grid[key]],
                mode="lines",
                name=label,
                line={"color": color, "width": 1.6},
            ),
            row=1,
            col=2,
        )

    fig.update_layout(template="plotly_white", legend_title="Fit")
    fig.update_xaxes(title_text=config.x_title, row=1, col=1)
    fig.update_yaxes(title_text=config.y_title, row=1, col=1)
    fig.update_xaxes(title_text=config.x_title, row=1, col=2)
    fig.update_yaxes(title_text="Delta chi_hat / chi_hat", row=1, col=2, rangemode="tozero")
    _write_figure_html(fig, path)


def save_dimensionless_fit_outputs(
    output_dir: Path,
    diag: list[dict[str, Any]],
    config: FitScaleConfig,
    *,
    max_rel_err: float | None,
    bootstrap_samples: int,
    seed: int,
) -> dict[str, Any]:
    fit_dir = output_dir / f"fits_{config.name}"
    flow_entries: dict[str, Any] = {}
    for flow_time in sorted({float(row["t_over_a2"]) for row in diag if config.x_key in row}):
        rows = [row for row in diag if np.isclose(float(row["t_over_a2"]), flow_time)]
        fit = build_dimensionless_fit_for_flow(
            rows,
            output_dir,
            config,
            max_rel_err=max_rel_err,
            bootstrap_samples=bootstrap_samples,
            seed=seed + int(round(1000.0 * flow_time)),
        )
        if fit is None:
            flow_entries[f"{flow_time:g}"] = {
                "status": "skipped",
                "reason": "fewer than three finite positive points after fit filters",
            }
            continue

        token = filename_token(f"{flow_time:g}")
        json_path = fit_dir / f"fit_t_over_a2_{token}.json"
        png_path = fit_dir / f"fit_t_over_a2_{token}.png"
        html_path = fit_dir / f"fit_t_over_a2_{token}.html"
        save_json(json_path, fit)
        fit_plot_paths = save_fit_png(png_path, fit, config)
        save_fit_html(html_path, fit, config)
        flow_entries[f"{flow_time:g}"] = {
            "status": "ok",
            "json": str(json_path),
            "png": str(fit_plot_paths.get("png", png_path)),
            "pdf": str(fit_plot_paths["pdf"]) if "pdf" in fit_plot_paths else None,
            "html": str(html_path),
            "plots": {fmt: str(fmt_path) for fmt, fmt_path in fit_plot_paths.items()},
            "n_points_fit": int(fit["n_points_fit"]),
            "n_points_available": int(fit["n_points_available"]),
        }

    summary = {
        "scale": config.name,
        "x_key": config.x_key,
        "y_key": config.y_key,
        "err_key": config.err_key,
        "fit_options": {
            "max_rel_err": float(max_rel_err) if max_rel_err is not None else None,
            "bootstrap_samples": int(bootstrap_samples),
            "seed": int(seed),
        },
        "flows": flow_entries,
    }
    summary_path = fit_dir / "summary.json"
    save_json(summary_path, summary)
    summary["summary_path"] = str(summary_path)
    return summary


def plot_html(
    path: Path,
    records: list[dict[str, Any]],
    *,
    x_key: str,
    y_key: str,
    err_key: str,
    title: str,
    x_title: str,
    y_title: str,
) -> None:
    go = _get_plotly()
    fig = go.Figure()
    for flow_time in sorted({float(record["t_over_a2"]) for record in records}):
        rows = sorted(
            [record for record in records if float(record["t_over_a2"]) == flow_time and y_key in record],
            key=lambda row: float(row[x_key]),
        )
        if not rows:
            continue
        fig.add_trace(
            go.Scatter(
                x=[float(row[x_key]) for row in rows],
                y=[float(row[y_key]) for row in rows],
                error_y={
                    "type": "data",
                    "array": [
                        float(row[err_key]) if row.get(err_key) is not None else 0.0
                        for row in rows
                    ],
                    "visible": True,
                },
                customdata=[
                    [float(row["R_mid"]), float(row["T_mid"]), float(row["t_over_a2"])]
                    for row in rows
                ],
                hovertemplate=(
                    "Rmid=%{customdata[0]:g}<br>"
                    "Tmid=%{customdata[1]:g}<br>"
                    "t_lat=%{customdata[2]:g}<br>"
                    f"{y_key}=%{{y:g}}<extra></extra>"
                ),
                mode="lines+markers",
                name=f"{flow_time:g}",
            )
        )

    fig.update_layout(
        title=title,
        template="plotly_white",
        xaxis_title=x_title,
        yaxis_title=y_title,
        legend_title="t_lat",
    )
    _write_figure_html(fig, path)


def plot_png(
    path: Path,
    records: list[dict[str, Any]],
    *,
    x_key: str,
    y_key: str,
    err_key: str,
    x_title: str,
    y_title: str,
) -> dict[str, Path]:
    fig, ax = plt.subplots(figsize=(4.8, 3.8), dpi=160)
    for flow_time in sorted({float(record["t_over_a2"]) for record in records}):
        rows = sorted(
            [record for record in records if float(record["t_over_a2"]) == flow_time and y_key in record],
            key=lambda row: float(row[x_key]),
        )
        if not rows:
            continue
        yerr = [
            float(row[err_key]) if row.get(err_key) is not None else 0.0
            for row in rows
        ]
        ax.errorbar(
            [float(row[x_key]) for row in rows],
            [float(row[y_key]) for row in rows],
            yerr=yerr,
            marker="o",
            markersize=2.5,
            linewidth=0.8,
            capsize=1.5,
            label=f"{flow_time:g}",
        )
    ax.set_xlabel(x_title)
    ax.set_ylabel(y_title)
    ax.legend(title=r"$t_{\mathrm{lat}}$", frameon=True, fontsize=8, title_fontsize=8)
    ax.tick_params(direction="in", top=True, right=True)
    fig.tight_layout()
    written = save_matplotlib_figure(fig, path)
    plt.close(fig)
    if not written:
        raise RuntimeError(f"No plot exports could be written for {path}.")
    return written


def save_plots(output_dir: Path, records: list[dict[str, Any]], *, t0_over_a2: float | None, r0_over_a: float | None) -> dict[str, str]:
    paths: dict[str, str] = {}
    diag = diagonal_records(records)
    if not diag:
        return paths

    raw_html = output_dir / "creutz_ratios_diagonal.html"
    raw_png = output_dir / "creutz_ratios_diagonal.png"
    plot_html(
        raw_html,
        diag,
        x_key="R_mid",
        y_key="chi",
        err_key="chi_err",
        title="Diagonal Creutz ratios from precomputed Wilson loops",
        x_title="R = T midpoint",
        y_title="chi_lattice",
    )
    raw_plot_paths = plot_png(raw_png, diag, x_key="R_mid", y_key="chi", err_key="chi_err", x_title="R = T midpoint", y_title="chi_lattice")
    paths["raw_diagonal_html"] = str(raw_html)
    for fmt, fmt_path in raw_plot_paths.items():
        paths[f"raw_diagonal_{fmt}"] = str(fmt_path)

    if t0_over_a2 is not None:
        dim_html = output_dir / "dimensionless_creutz_ratios_t0_diagonal.html"
        dim_png = output_dir / "dimensionless_creutz_ratios_t0_diagonal.png"
        plot_html(
            dim_html,
            diag,
            x_key="r_hat_t0",
            y_key="chi_hat_t0",
            err_key="chi_hat_t0_err",
            title="Dimensionless diagonal Creutz ratios",
            x_title="r / sqrt(8 t0)",
            y_title="chi_hat = chi_lattice * 8 t0/a^2",
        )
        dim_plot_paths = plot_png(
            dim_png,
            diag,
            x_key="r_hat_t0",
            y_key="chi_hat_t0",
            err_key="chi_hat_t0_err",
            x_title=r"$\hat r = r/\sqrt{8t_0}$",
            y_title=r"$\hat\chi = \chi_{\rm lat}\,8t_0/a^2$",
        )
        paths["dimensionless_t0_diagonal_html"] = str(dim_html)
        for fmt, fmt_path in dim_plot_paths.items():
            paths[f"dimensionless_t0_diagonal_{fmt}"] = str(fmt_path)

    if r0_over_a is not None:
        dim_html = output_dir / "dimensionless_creutz_ratios_r0_diagonal.html"
        dim_png = output_dir / "dimensionless_creutz_ratios_r0_diagonal.png"
        plot_html(
            dim_html,
            diag,
            x_key="r_hat_r0",
            y_key="chi_hat_r0",
            err_key="chi_hat_r0_err",
            title="r0-scaled diagonal Creutz ratios",
            x_title="r / r0",
            y_title="chi_hat = chi_lattice * (r0/a)^2",
        )
        dim_plot_paths = plot_png(
            dim_png,
            diag,
            x_key="r_hat_r0",
            y_key="chi_hat_r0",
            err_key="chi_hat_r0_err",
            x_title=r"$\hat r = r/r_0$",
            y_title=r"$\hat\chi = \chi_{\rm lat}(r_0/a)^2$",
        )
        paths["dimensionless_r0_diagonal_html"] = str(dim_html)
        for fmt, fmt_path in dim_plot_paths.items():
            paths[f"dimensionless_r0_diagonal_{fmt}"] = str(fmt_path)

    return paths


def output_dir_for(analysis_dir: Path, output_root: Path | None, n_groups: int) -> Path:
    if output_root is None:
        return analysis_dir / "creutz_ratio_analysis"
    root = output_root.expanduser().resolve()
    if n_groups == 1:
        return root
    return root / analysis_dir.name


def analyze_one(
    analysis_dir: Path,
    output_dir: Path,
    *,
    t0_over_a2: float | None,
    r0_over_a: float | None,
    use_bootstrap: bool,
    make_fits: bool,
    fit_max_rel_err: float | None,
    fit_bootstrap_samples: int,
    fit_seed: int,
) -> dict[str, Any]:
    rows = load_wilson_loop_stats(analysis_dir)
    records = compute_creutz_records(
        analysis_dir,
        output_dir,
        rows,
        t0_over_a2=t0_over_a2,
        r0_over_a=r0_over_a,
        use_bootstrap=use_bootstrap,
        diagonal_only=False,
    )
    diag = diagonal_records(records)

    output_dir.mkdir(parents=True, exist_ok=True)
    table_path = output_dir / "creutz_ratios.dat"
    diagonal_table_path = output_dir / "creutz_ratios_diagonal.dat"
    write_table(table_path, records)
    write_table(diagonal_table_path, diag)
    plot_paths = save_plots(output_dir, records, t0_over_a2=t0_over_a2, r0_over_a=r0_over_a)
    fit_summaries: dict[str, Any] = {}
    if make_fits:
        if t0_over_a2 is not None:
            scale_factor = float(8.0 * t0_over_a2)
            fit_summaries["t0"] = save_dimensionless_fit_outputs(
                output_dir,
                diag,
                FitScaleConfig(
                    name="t0",
                    x_key="r_hat_t0",
                    y_key="chi_hat_t0",
                    err_key="chi_hat_t0_err",
                    scale_factor=scale_factor,
                    x_denominator=math.sqrt(scale_factor),
                    x_title="r / sqrt(8 t0)",
                    y_title="chi_hat = chi_lattice * 8 t0/a^2",
                    x_label=r"$\hat r = r/\sqrt{8t_0}$",
                    y_label=r"$\hat\chi$",
                ),
                max_rel_err=fit_max_rel_err,
                bootstrap_samples=fit_bootstrap_samples,
                seed=fit_seed,
            )
        if r0_over_a is not None:
            scale_factor = float(r0_over_a) ** 2
            fit_summaries["r0"] = save_dimensionless_fit_outputs(
                output_dir,
                diag,
                FitScaleConfig(
                    name="r0",
                    x_key="r_hat_r0",
                    y_key="chi_hat_r0",
                    err_key="chi_hat_r0_err",
                    scale_factor=scale_factor,
                    x_denominator=float(r0_over_a),
                    x_title="r / r0",
                    y_title="chi_hat = chi_lattice * (r0/a)^2",
                    x_label=r"$\hat r = r/r_0$",
                    y_label=r"$\hat\chi$",
                ),
                max_rel_err=fit_max_rel_err,
                bootstrap_samples=fit_bootstrap_samples,
                seed=fit_seed,
            )

    source_summary_path = analysis_dir / "summary.json"
    source_summary = None
    if source_summary_path.exists():
        with source_summary_path.open("r", encoding="utf-8") as handle:
            source_summary = json.load(handle)

    summary = {
        "schema_version": 1,
        "source_analysis_dir": str(analysis_dir),
        "source_wilson_loop_table": str(analysis_dir / STATS_FILENAME),
        "source_summary": source_summary,
        "t0_over_a2": t0_over_a2,
        "r0_over_a": r0_over_a,
        "used_wilson_bootstrap_samples": bool(use_bootstrap),
        "n_wilson_loop_stats": len(rows),
        "n_creutz_ratios": len(records),
        "n_diagonal_creutz_ratios": len(diag),
        "flow_times": sorted({record["t_over_a2"] for record in records}),
        "tables": {
            "all": str(table_path),
            "diagonal": str(diagonal_table_path),
        },
        "plots": plot_paths,
        "fits": fit_summaries,
        "records": records,
        "saved_at": datetime.now(timezone.utc).isoformat(),
    }
    summary_path = output_dir / "summary.json"
    save_json(summary_path, summary)
    summary["summary_path"] = str(summary_path)
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compute Creutz ratios from precomputed gradient_flow_wtemp_analysis "
            "Wilson-loop means and bootstrap samples."
        )
    )
    parser.add_argument(
        "paths",
        nargs="+",
        type=Path,
        help=(
            "One or more gradient_flow_wtemp_analysis directories, wilson_loop_stats.dat "
            "files, or roots containing such directories."
        ),
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        help=(
            "Output directory. For multiple discovered analyses, one subdirectory per "
            "input analysis is created. Default: <analysis>/creutz_ratio_analysis."
        ),
    )
    parser.add_argument(
        "--t0-over-a2",
        type=float,
        help=(
            "Reference flow time in lattice units. Enables paper-style dimensionless "
            "columns and plots: chi_hat=chi_lattice*8*t0/a^2 and r_hat=r/sqrt(8*t0)."
        ),
    )
    parser.add_argument(
        "--r0-over-a",
        type=float,
        help=(
            "Sommer scale in lattice units. Enables r0-scaled columns and plots: "
            "chi_hat=chi_lattice*(r0/a)^2 and r_hat=r/r0."
        ),
    )
    parser.add_argument(
        "--no-bootstrap",
        action="store_true",
        help="Do not load Wilson-loop bootstrap arrays; use independent-error propagation instead.",
    )
    parser.add_argument(
        "--no-fits",
        action="store_true",
        help="Skip per-flow dimensionless fit plots.",
    )
    parser.add_argument(
        "--fit-max-rel-err",
        type=float,
        help="Maximum relative point error included in interpolations. Default: no relative-error cut.",
    )
    parser.add_argument(
        "--fit-bootstrap-samples",
        type=int,
        default=200,
        help="Maximum bootstrap replicas used for the fit statistical band. Default: 200.",
    )
    parser.add_argument(
        "--fit-seed",
        type=int,
        default=12345,
        help="Random seed for Gaussian fit replicas when bootstrap arrays are unavailable.",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List discovered gradient_flow_wtemp_analysis directories and exit.",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    if args.t0_over_a2 is not None and args.t0_over_a2 <= 0:
        parser.error("--t0-over-a2 must be positive")
    if args.r0_over_a is not None and args.r0_over_a <= 0:
        parser.error("--r0-over-a must be positive")
    if args.fit_max_rel_err is not None and args.fit_max_rel_err <= 0:
        parser.error("--fit-max-rel-err must be positive")
    if args.fit_bootstrap_samples < 1:
        parser.error("--fit-bootstrap-samples must be >= 1")

    try:
        analysis_dirs = discover_analysis_dirs(args.paths)
    except (FileNotFoundError, ValueError) as exc:
        parser.error(str(exc))

    if not analysis_dirs:
        parser.error("No gradient_flow_wtemp_analysis directories with wilson_loop_stats.dat were found.")

    if args.list:
        for analysis_dir in analysis_dirs:
            print(analysis_dir)
        return 0

    for index, analysis_dir in enumerate(analysis_dirs, start=1):
        output_dir = output_dir_for(analysis_dir, args.output_dir, len(analysis_dirs))
        if len(analysis_dirs) > 1:
            print(f"[{index}/{len(analysis_dirs)}] {analysis_dir}")
        summary = analyze_one(
            analysis_dir,
            output_dir,
            t0_over_a2=args.t0_over_a2,
            r0_over_a=args.r0_over_a,
            use_bootstrap=not args.no_bootstrap,
            make_fits=not args.no_fits and (args.t0_over_a2 is not None or args.r0_over_a is not None),
            fit_max_rel_err=args.fit_max_rel_err,
            fit_bootstrap_samples=args.fit_bootstrap_samples,
            fit_seed=args.fit_seed,
        )
        print(
            f"Saved {summary['n_creutz_ratios']} Creutz ratios "
            f"({summary['n_diagonal_creutz_ratios']} diagonal) to {output_dir}"
        )
        for label, path in sorted(summary["plots"].items()):
            print(f"  {label}: {path}")
        for scale_name, fit_summary in sorted(summary["fits"].items()):
            print(f"  {scale_name} fits: {fit_summary['summary_path']}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
