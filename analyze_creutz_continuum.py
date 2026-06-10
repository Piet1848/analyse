#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import os
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from analyze_precomputed_creutz_ratios import (
    FitScaleConfig,
    filtered_fit_rows,
    fit_model_curves,
    load_dimensionless_bootstrap_matrix,
)
from finalized_analysis_helpers import _get_plotly, _write_figure_html, save_matplotlib_figure


SUMMARY_FILENAME = "summary.json"
DEFAULT_OUTPUT_ROOT = Path("../data/creutz_continuum_analysis")
FIT_TERM_NAMES = ("c00", "c20", "c40", "c01", "c21", "c02")
SEPARATE_LINEAR_TERM_NAMES = ("c0", "c2")
SEPARATE_QUADRATIC_TERM_NAMES = ("c0", "c2", "c4")
FIT_MODES = ("separate-linear", "separate-quadratic", "combined-risch")
DEFAULT_INPUT_ROOT = Path("../data/gradient_flow_wtemp_analysis")
DEFAULT_BETAS = (2.3, 2.4, 2.5, 2.6)
DEFAULT_EPS1 = (0.0,)
DEFAULT_FLOW_TIMES = (0.0, 0.25, 0.5, 0.75, 1.0)
DEFAULT_T0_OVER_A2_BY_BETA = {
    2.3: 0.5689746602537821,
    2.4: 0.9681815767085601,
    2.5: 1.8481372605040205,
    2.6: 3.565414260750965,
}


@dataclass(frozen=True)
class ContinuumInput:
    analysis_dir: Path
    summary: dict[str, Any]
    beta: float | None
    inferred_epsilon: float | None
    scale_over_a: float | None
    summary_t0_over_a2: float | None
    fit_summary_path: Path
    flow_paths: dict[float, Path]


@dataclass(frozen=True)
class PredictionBand:
    y: np.ndarray
    y_low: np.ndarray
    y_high: np.ndarray
    source: str


def filename_token(value: float | str) -> str:
    return str(value).replace("-", "m").replace(".", "p").replace(",", "_")


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


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def finite_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if np.isfinite(out) else None


def row_t_lat(row: dict[str, Any]) -> float | None:
    value = finite_float(row.get("t_lat"))
    if value is not None:
        return value
    value = finite_float(row.get("t_over_a2"))
    if value is not None:
        return value
    eight_t = finite_float(row.get("eight_t_over_a2"))
    return None if eight_t is None else eight_t / 8.0


def t_lat_to_eight_t_over_a2(t_lat: float) -> float:
    return 8.0 * float(t_lat)


def t_lat_label(t_lat: float) -> str:
    return f"t_lat={float(t_lat):g}"


def t_lat_legend_label(t_lat: float, suffix: str | None = None) -> str:
    label = rf"t_{{\mathrm{{lat}}}}={float(t_lat):g}"
    if suffix:
        label += rf"\ \mathrm{{{suffix}}}"
    return f"${label}$"


def continuum_t_lat_legend_label(t_lat: float) -> str:
    return rf"$a=0,\ t_{{\mathrm{{lat}}}}={float(t_lat):g}$"


def prompt_float(prompt: str, default: float | None = None) -> float:
    suffix = f" [{default:g}]" if default is not None and np.isfinite(default) else ""
    while True:
        raw = input(f"{prompt}{suffix}: ").strip()
        if not raw and default is not None:
            return float(default)
        try:
            return float(raw)
        except ValueError:
            print(f"Could not parse a number from {raw!r}.")


def discover_analysis_dirs(paths: list[Path]) -> list[Path]:
    seen: set[Path] = set()
    result: list[Path] = []
    for raw_path in paths:
        path = raw_path.expanduser().resolve()
        candidates: list[Path] = []
        if path.is_file() and path.name == SUMMARY_FILENAME:
            candidates = [path.parent]
        elif (path / SUMMARY_FILENAME).exists() and any((path / name / SUMMARY_FILENAME).exists() for name in ("fits_r0", "fits_t0")):
            candidates = [path]
        elif path.is_dir():
            candidates = sorted(
                item.parent
                for item in path.rglob(SUMMARY_FILENAME)
                if item.parent.name == "creutz_ratio_analysis"
            )
        else:
            raise FileNotFoundError(f"No such input path: {path}")

        for candidate in candidates:
            if candidate not in seen:
                seen.add(candidate)
                result.append(candidate)
    return result


def group_signature_value(summary: dict[str, Any], name: str) -> Any:
    source_summary = summary.get("source_summary") or {}
    for row in source_summary.get("group_signature", []) or []:
        if row.get("name") == name:
            return row.get("value")
    return None


def infer_from_path(path: Path, key: str) -> float | None:
    patterns = {
        "beta": r"beta_([0-9pm]+)",
        "epsilon1": r"epsilon1_([0-9pm]+)",
    }
    match = re.search(patterns[key], str(path))
    if not match:
        return None
    token = match.group(1).replace("p", ".").replace("m", "-")
    return finite_float(token)


def fit_flow_paths(fit_summary: dict[str, Any], fit_summary_path: Path) -> dict[float, Path]:
    paths: dict[float, Path] = {}
    for key, row in (fit_summary.get("flows") or {}).items():
        if not isinstance(row, dict) or row.get("status") != "ok" or not row.get("json"):
            continue
        flow_time = finite_float(key)
        if flow_time is None:
            continue
        path = Path(row["json"])
        if not path.is_absolute():
            path = fit_summary_path.parent / path
        paths[flow_time] = path
    return paths


def load_continuum_input(analysis_dir: Path, scale: str) -> ContinuumInput | None:
    summary_path = analysis_dir / SUMMARY_FILENAME
    fit_summary_path = analysis_dir / f"fits_{scale}" / SUMMARY_FILENAME
    if not summary_path.exists() or not fit_summary_path.exists():
        return None

    summary = load_json(summary_path)
    fit_summary = load_json(fit_summary_path)
    flow_paths = fit_flow_paths(fit_summary, fit_summary_path)
    if not flow_paths:
        return None

    beta = finite_float(group_signature_value(summary, "MetropolisParams.beta"))
    epsilon = finite_float(group_signature_value(summary, "MetropolisParams.epsilon1"))
    beta = beta if beta is not None else infer_from_path(analysis_dir, "beta")
    epsilon = epsilon if epsilon is not None else infer_from_path(analysis_dir, "epsilon1")

    summary_t0_over_a2 = finite_float(summary.get("t0_over_a2"))
    scale_over_a: float | None
    if scale == "r0":
        scale_over_a = finite_float(summary.get("r0_over_a"))
        if scale_over_a is None:
            return None
    else:
        scale_over_a = math.sqrt(8.0 * summary_t0_over_a2) if summary_t0_over_a2 is not None else None

    return ContinuumInput(
        analysis_dir=analysis_dir,
        summary=summary,
        beta=beta,
        inferred_epsilon=epsilon,
        scale_over_a=scale_over_a,
        summary_t0_over_a2=summary_t0_over_a2,
        fit_summary_path=fit_summary_path,
        flow_paths=flow_paths,
    )


def nearest_flow_time(available: list[float], requested: float, *, atol: float = 1e-9) -> float | None:
    if not available:
        return None
    best = min(available, key=lambda value: abs(value - requested))
    return best if abs(best - requested) <= atol else None


def common_flow_times(inputs: list[ContinuumInput]) -> list[float]:
    if not inputs:
        return []
    common = set(inputs[0].flow_paths)
    for item in inputs[1:]:
        common &= set(item.flow_paths)
    return sorted(common)


def float_in_list(value: float | None, allowed: list[float], *, atol: float = 1e-8) -> bool:
    if value is None:
        return False
    return any(abs(float(value) - float(item)) <= atol for item in allowed)


def filter_inputs_by_beta(inputs: list[ContinuumInput], betas: list[float]) -> list[ContinuumInput]:
    return [item for item in inputs if float_in_list(item.beta, betas)]


def fit_exclusion_reason(row: dict[str, Any], excluded_betas: list[float]) -> str | None:
    if excluded_betas and float_in_list(finite_float(row.get("beta")), excluded_betas):
        return "beta excluded from continuum fit"
    return None


def mark_fit_usage(rows: list[dict[str, Any]], excluded_betas: list[float]) -> list[dict[str, Any]]:
    fit_rows: list[dict[str, Any]] = []
    for row in rows:
        reason = fit_exclusion_reason(row, excluded_betas)
        row["used_in_fit"] = reason is None
        row["fit_exclusion_reason"] = reason
        if reason is None:
            fit_rows.append(row)
    return fit_rows


def resolve_flow_times(args: argparse.Namespace, inputs: list[ContinuumInput]) -> list[float]:
    common = common_flow_times(inputs)
    if args.eight_t_over_a2 is not None:
        requested = [float(value) / 8.0 for value in args.eight_t_over_a2]
    elif args.flow_time is not None:
        requested = [float(value) for value in args.flow_time]
    else:
        requested = list(DEFAULT_FLOW_TIMES)

    if not requested:
        raise ValueError("No common t_lat exists across all selected inputs.")

    selected: list[float] = []
    for value in requested:
        flow = nearest_flow_time(common, value)
        if flow is None:
            available_text = ", ".join(f"{item:g}" for item in common)
            raise ValueError(f"Requested t_lat={value:g} is not common to all inputs. Common values: {available_text}")
        if flow not in selected:
            selected.append(flow)
    return selected


def resolve_r_hat(args: argparse.Namespace) -> float:
    if args.r_hat is not None:
        return float(args.r_hat)
    if os.isatty(0):
        return prompt_float("Target r_hat")
    raise ValueError("--r-hat is required in non-interactive runs")


def parse_value_map(tokens: list[str] | None, option_name: str) -> dict[float, float]:
    result: dict[float, float] = {}
    for token in tokens or []:
        if "=" not in token:
            raise ValueError(f"{option_name} entries must use BETA=VALUE syntax: {token!r}")
        raw_beta, raw_value = token.split("=", 1)
        beta = finite_float(raw_beta)
        value = finite_float(raw_value)
        if beta is None or value is None:
            raise ValueError(f"{option_name} entry contains a non-finite number: {token!r}")
        result[beta] = value
    return result


def resolve_eps_map_from_args(args: argparse.Namespace, betas: list[float]) -> dict[float, float]:
    if args.eps_map:
        return parse_value_map(args.eps_map, "--eps-map")
    eps_values = [float(value) for value in args.eps1]
    if len(eps_values) == 1:
        return {float(beta): eps_values[0] for beta in betas}
    if len(eps_values) == len(betas):
        return {float(beta): float(eps) for beta, eps in zip(betas, eps_values, strict=True)}
    raise ValueError(
        "--eps1 must contain either one value applied to all selected betas, "
        "or the same number of values as --betas."
    )


def resolve_t0_map_from_args(args: argparse.Namespace, betas: list[float]) -> dict[float, float]:
    result: dict[float, float] = {}
    for beta in betas:
        default = lookup_by_beta(DEFAULT_T0_OVER_A2_BY_BETA, beta)
        if default is not None:
            result[float(beta)] = float(default)
    result.update(parse_value_map(args.t0_map, "--t0-map"))
    return result


def lookup_by_beta(mapping: dict[float, float], beta: float | None, *, atol: float = 1e-9) -> float | None:
    if beta is None:
        return None
    for key, value in mapping.items():
        if abs(key - beta) <= atol:
            return value
    return None


def beta_sort_key(value: float | None) -> tuple[float, str]:
    if value is None:
        return math.inf, "unknown"
    return value, f"{value:g}"


def group_inputs_by_beta(inputs: list[ContinuumInput]) -> dict[float | None, list[ContinuumInput]]:
    groups: dict[float | None, list[ContinuumInput]] = {}
    for item in inputs:
        groups.setdefault(item.beta, []).append(item)
    return groups


def resolve_fixed_eps_selection(
    inputs: list[ContinuumInput],
    eps_map: dict[float, float],
    *,
    prompt: bool,
    eps_tol: float,
) -> tuple[list[ContinuumInput], dict[Path, float | None], dict[str, Any], list[dict[str, Any]]]:
    selected: list[ContinuumInput] = []
    eps_by_dir: dict[Path, float | None] = {}
    chosen_by_beta: dict[str, float | None] = {}
    skipped: list[dict[str, Any]] = []

    groups = group_inputs_by_beta(inputs)
    for beta, group in sorted(groups.items(), key=lambda pair: beta_sort_key(pair[0])):
        available = sorted({item.inferred_epsilon for item in group if item.inferred_epsilon is not None})
        chosen = lookup_by_beta(eps_map, beta)
        beta_label = "unknown" if beta is None else f"{beta:g}"
        if chosen is None:
            if len(available) == 1:
                chosen = available[0]
            elif prompt and os.isatty(0) and available:
                print(f"Available epsilon1 values for beta={beta_label}: " + ", ".join(f"{value:g}" for value in available))
                chosen = prompt_float(f"Fixed-trajectory epsilon for beta={beta_label}", default=available[0])
            elif not available:
                chosen = None
            else:
                values = ", ".join(f"{value:g}" for value in available)
                raise ValueError(
                    f"Multiple epsilon1 values were found for beta={beta_label}: {values}. "
                    "Pass --eps-map BETA=EPS to select the fixed-epsilon trajectory."
                )

        chosen_by_beta[beta_label] = chosen
        for item in group:
            if chosen is not None and item.inferred_epsilon is not None and abs(item.inferred_epsilon - chosen) > eps_tol:
                skipped.append(
                    {
                        "analysis_dir": str(item.analysis_dir),
                        "beta": item.beta,
                        "epsilon_inferred": item.inferred_epsilon,
                        "epsilon_fixed": chosen,
                        "reason": "epsilon does not match selected fixed-epsilon trajectory",
                    }
                )
                continue
            selected.append(item)
            eps_by_dir[item.analysis_dir] = chosen

    return selected, eps_by_dir, {"by_beta": chosen_by_beta, "tolerance": eps_tol}, skipped


def resolve_t0_values(
    inputs: list[ContinuumInput],
    t0_map: dict[float, float],
    *,
    prompt: bool,
) -> tuple[dict[Path, float], dict[str, float]]:
    values_by_dir: dict[Path, float] = {}
    chosen_by_beta: dict[str, float] = {}
    groups = group_inputs_by_beta(inputs)

    for beta, group in sorted(groups.items(), key=lambda pair: beta_sort_key(pair[0])):
        beta_label = "unknown" if beta is None else f"{beta:g}"
        mapped = lookup_by_beta(t0_map, beta)
        defaults = sorted({item.summary_t0_over_a2 for item in group if item.summary_t0_over_a2 is not None})
        chosen = mapped if mapped is not None else (defaults[0] if defaults else None)

        if chosen is None and prompt and os.isatty(0):
            chosen = prompt_float(f"t0/a^2 for beta={beta_label}")
        if chosen is None:
            raise ValueError(f"Missing t0/a^2 for beta={beta_label}. Pass --t0-map {beta_label}=T0_OVER_A2.")
        if chosen <= 0:
            raise ValueError(f"t0/a^2 must be positive for beta={beta_label}; got {chosen:g}.")

        chosen_by_beta[beta_label] = float(chosen)
        for item in group:
            values_by_dir[item.analysis_dir] = float(chosen)

    return values_by_dir, chosen_by_beta


def interpolate_fit_at_r_hat(fit_json_path: Path, r_hat: float) -> dict[str, float] | None:
    fit = load_json(fit_json_path)
    grid = fit.get("grid") or {}
    xs = np.asarray(grid.get("x", []), dtype=float)
    finite = np.isfinite(xs)
    if xs.size == 0 or not np.any(finite):
        return None
    x_min = float(np.nanmin(xs))
    x_max = float(np.nanmax(xs))
    if r_hat < x_min or r_hat > x_max:
        return None

    def interp(key: str) -> float | None:
        values = np.asarray(grid.get(key, []), dtype=float)
        if values.size != xs.size:
            return None
        valid = np.isfinite(xs) & np.isfinite(values)
        if np.count_nonzero(valid) < 2:
            return None
        order = np.argsort(xs[valid])
        return float(np.interp(r_hat, xs[valid][order], values[valid][order]))

    stat_err = interp("stat_err")
    sys_err = interp("sys_err")
    total_err = interp("total_err")
    if stat_err is not None and sys_err is not None:
        total_err = float(math.hypot(stat_err, sys_err))

    return {
        "chi_hat": interp("average"),
        "stat_err": stat_err,
        "sys_err": sys_err,
        "total_err": total_err,
        "r_hat_min": x_min,
        "r_hat_max": x_max,
    }


def sigma_from_row(row: dict[str, Any]) -> float:
    total = finite_float(row.get("total_err"))
    if total is not None and total > 0:
        return total
    stat = finite_float(row.get("stat_err")) or 0.0
    sys = finite_float(row.get("sys_err")) or 0.0
    total = math.hypot(stat, sys)
    return total if total > 0 and np.isfinite(total) else 1.0


def build_rows(
    inputs: list[ContinuumInput],
    *,
    scale: str,
    r_hat: float,
    flow_times: list[float],
    fixed_eps_by_dir: dict[Path, float | None],
    t0_by_dir: dict[Path, float],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    for item in inputs:
        t0_over_a2 = t0_by_dir[item.analysis_dir]
        if item.scale_over_a is None or item.scale_over_a <= 0:
            skipped.append({"analysis_dir": str(item.analysis_dir), "reason": f"missing positive {scale}/a scale"})
            continue
        ahat2 = 1.0 / (float(item.scale_over_a) ** 2)
        for flow_time in flow_times:
            path = item.flow_paths.get(flow_time)
            if path is None:
                skipped.append({"analysis_dir": str(item.analysis_dir), "t_over_a2": flow_time, "reason": "missing flow fit"})
                continue
            value = interpolate_fit_at_r_hat(path, r_hat)
            if value is None or value["chi_hat"] is None:
                skipped.append(
                    {
                        "analysis_dir": str(item.analysis_dir),
                        "t_over_a2": flow_time,
                        "reason": "target r_hat outside fit grid",
                        "fit_json": str(path),
                    }
                )
                continue

            smearing_strength = t_lat_to_eight_t_over_a2(flow_time)
            tau = flow_time / t0_over_a2
            rows.append(
                {
                    "analysis_dir": str(item.analysis_dir),
                    "fit_json": str(path),
                    "beta": item.beta,
                    "epsilon_inferred": item.inferred_epsilon,
                    "epsilon_fixed": fixed_eps_by_dir.get(item.analysis_dir),
                    "scale_over_a": item.scale_over_a,
                    "t0_over_a2": float(t0_over_a2),
                    "ahat2": float(ahat2),
                    "r_hat": float(r_hat),
                    "t_over_a2": float(flow_time),
                    "t_lat": float(flow_time),
                    "eight_t_over_a2": float(smearing_strength),
                    "tau": float(tau),
                    **value,
                }
            )

    return sorted(rows, key=lambda row: (float(row_t_lat(row) or 0.0), float(row["ahat2"]))), skipped


def bootstrap_scale_config(item: ContinuumInput, scale: str, t0_over_a2: float) -> FitScaleConfig | None:
    if scale == "r0":
        if item.scale_over_a is None:
            return None
        scale_factor = float(item.scale_over_a) ** 2
        return FitScaleConfig(
            name="r0",
            x_key="r_hat_r0",
            y_key="chi_hat_r0",
            err_key="chi_hat_r0_err",
            scale_factor=scale_factor,
            x_denominator=float(item.scale_over_a),
            x_title="r / r0",
            y_title="chi_hat = chi_lattice * (r0/a)^2",
            x_label=r"$\hat r = r/r_0$",
            y_label=r"$\hat\chi$",
        )

    scale_factor = 8.0 * float(t0_over_a2)
    return FitScaleConfig(
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
    )


def diagonal_records_for_flow(item: ContinuumInput, flow_time: float, config: FitScaleConfig) -> list[dict[str, Any]]:
    records = item.summary.get("records") or []
    selected: list[dict[str, Any]] = []
    for row in records:
        if config.x_key not in row or config.y_key not in row or config.err_key not in row:
            continue
        if not np.isclose(float(row.get("t_over_a2", float("nan"))), flow_time):
            continue
        if not np.isclose(float(row.get("R_mid", float("nan"))), float(row.get("T_mid", float("nan")))):
            continue
        selected.append(row)
    return sorted(selected, key=lambda row: float(row[config.x_key]))


def interpolated_bootstrap_values(
    item: ContinuumInput,
    *,
    scale: str,
    flow_time: float,
    r_hat: float,
    fit_json_path: Path,
    t0_over_a2: float,
    max_samples: int,
    seed: int,
) -> tuple[np.ndarray | None, dict[str, Any]]:
    config = bootstrap_scale_config(item, scale, t0_over_a2)
    if config is None:
        return None, {"source": "unavailable", "reason": "missing scale information"}

    fit_json = load_json(fit_json_path)
    max_rel_err = finite_float((fit_json.get("fit_options") or {}).get("max_rel_err"))
    flow_rows = diagonal_records_for_flow(item, flow_time, config)
    fit_rows = filtered_fit_rows(
        flow_rows,
        x_key=config.x_key,
        y_key=config.y_key,
        err_key=config.err_key,
        max_rel_err=max_rel_err,
    )
    if len(fit_rows) < 3:
        return None, {"source": "unavailable", "reason": "fewer than three filtered interpolation nodes"}

    x = np.asarray([float(row[config.x_key]) for row in fit_rows], dtype=float)
    y = np.asarray([float(row[config.y_key]) for row in fit_rows], dtype=float)
    sigma = np.asarray([float(row[config.err_key]) for row in fit_rows], dtype=float)
    bootstrap_matrix = load_dimensionless_bootstrap_matrix(
        fit_rows,
        item.analysis_dir,
        scale_factor=config.scale_factor,
    )
    source = "creutz_bootstrap"
    if bootstrap_matrix is None:
        rng = np.random.default_rng(seed)
        bootstrap_matrix = rng.normal(y[:, None], sigma[:, None], size=(x.size, max(1, int(max_samples))))
        source = "gaussian_from_creutz_errors"

    n_available = int(bootstrap_matrix.shape[1])
    n_samples = min(max(1, int(max_samples)), n_available)
    if n_samples < n_available:
        indices = np.linspace(0, n_available - 1, n_samples, dtype=int)
    else:
        indices = np.arange(n_available, dtype=int)

    values = np.full(indices.size, np.nan, dtype=float)
    target = np.asarray([float(r_hat)], dtype=float)
    for out_index, sample_index in enumerate(indices):
        y_replica = np.asarray(bootstrap_matrix[:, sample_index], dtype=float)
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
            fitted = fit_model_curves(x[mask], y_replica[mask], sigma[mask], target)
        except Exception:
            continue
        average = np.asarray(fitted["average"], dtype=float)
        if average.size == 1 and np.isfinite(average[0]):
            values[out_index] = float(average[0])

    return values, {
        "source": source,
        "samples_available": n_available,
        "samples_requested": int(max_samples),
        "samples_returned": int(values.size),
        "finite_samples": int(np.count_nonzero(np.isfinite(values))),
        "n_interpolation_nodes": int(len(fit_rows)),
    }


def build_interpolated_bootstrap_matrix(
    rows: list[dict[str, Any]],
    inputs_by_dir: dict[str, ContinuumInput],
    *,
    scale: str,
    max_samples: int,
    seed: int,
) -> tuple[np.ndarray | None, list[dict[str, Any]]]:
    samples_by_row: list[np.ndarray] = []
    metadata: list[dict[str, Any]] = []
    for row_index, row in enumerate(rows):
        analysis_dir = str(row["analysis_dir"])
        item = inputs_by_dir.get(analysis_dir)
        if item is None:
            metadata.append({"row_index": row_index, "source": "unavailable", "reason": "missing input object"})
            return None, metadata

        values, info = interpolated_bootstrap_values(
            item,
            scale=scale,
            flow_time=float(row["t_over_a2"]),
            r_hat=float(row["r_hat"]),
            fit_json_path=Path(str(row["fit_json"])),
            t0_over_a2=float(row["t0_over_a2"]),
            max_samples=max_samples,
            seed=seed + row_index,
        )
        info["row_index"] = row_index
        info["analysis_dir"] = analysis_dir
        info["t_over_a2"] = float(row["t_over_a2"])
        metadata.append(info)
        if values is None:
            return None, metadata
        samples_by_row.append(values)

    if not samples_by_row:
        return None, metadata
    n_samples = min(sample.size for sample in samples_by_row)
    matrix = np.vstack([sample[:n_samples] for sample in samples_by_row])
    finite_columns = np.all(np.isfinite(matrix), axis=0)
    matrix = matrix[:, finite_columns]
    if matrix.shape[1] < 2:
        return None, metadata
    return matrix, metadata


def combined_design(ahat2: np.ndarray, tau: np.ndarray) -> np.ndarray:
    return np.column_stack(
        [
            np.ones_like(ahat2),
            ahat2,
            ahat2**2,
            tau,
            ahat2 * tau,
            tau**2,
        ]
    )


def separate_design(ahat2: np.ndarray, order: int) -> np.ndarray:
    if order == 1:
        return np.column_stack([np.ones_like(ahat2), ahat2])
    if order == 2:
        return np.column_stack([np.ones_like(ahat2), ahat2, ahat2**2])
    raise ValueError("order must be 1 or 2")


def separate_term_names(order: int) -> tuple[str, ...]:
    if order == 1:
        return SEPARATE_LINEAR_TERM_NAMES
    if order == 2:
        return SEPARATE_QUADRATIC_TERM_NAMES
    raise ValueError("order must be 1 or 2")


def fit_mode_order(mode: str) -> int:
    if mode == "separate-linear":
        return 1
    if mode == "separate-quadratic":
        return 2
    raise ValueError(f"{mode!r} is not a separate-smearing fit mode")


def design_for_rows(rows: list[dict[str, Any]]) -> np.ndarray:
    ahat2 = np.asarray([float(row["ahat2"]) for row in rows], dtype=float)
    tau = np.asarray([float(row["tau"]) for row in rows], dtype=float)
    return combined_design(ahat2, tau)


def regularized_covariance(matrix: np.ndarray) -> np.ndarray:
    covariance = 0.5 * (np.asarray(matrix, dtype=float) + np.asarray(matrix, dtype=float).T)
    diagonal = np.diag(covariance)
    finite_positive = diagonal[np.isfinite(diagonal) & (diagonal > 0)]
    scale = float(np.max(finite_positive)) if finite_positive.size else 1.0
    jitter = max(scale * 1e-12, 1e-15)
    return covariance + np.eye(covariance.shape[0]) * jitter


def block_statistical_covariance(rows: list[dict[str, Any]], bootstrap_y: np.ndarray) -> np.ndarray:
    n_rows = len(rows)
    covariance = np.zeros((n_rows, n_rows), dtype=float)
    groups: dict[str, list[int]] = {}
    for index, row in enumerate(rows):
        groups.setdefault(str(row["analysis_dir"]), []).append(index)

    for indices in groups.values():
        block = np.asarray(bootstrap_y[indices, :], dtype=float)
        if block.shape[1] < 2:
            continue
        block_cov = np.cov(block, ddof=1)
        if len(indices) == 1:
            covariance[indices[0], indices[0]] = float(block_cov)
        else:
            covariance[np.ix_(indices, indices)] = np.asarray(block_cov, dtype=float)
    return covariance


def solve_linear_ansatz(
    design: np.ndarray,
    y: np.ndarray,
    *,
    sigma: np.ndarray | None = None,
    covariance_y: np.ndarray | None = None,
) -> dict[str, Any]:
    if covariance_y is not None:
        cov_y = regularized_covariance(covariance_y)
        precision = np.linalg.pinv(cov_y)
        normal = design.T @ precision @ design
        covariance = np.linalg.pinv(normal)
        coeffs = covariance @ design.T @ precision @ y
        residuals = y - design @ coeffs
        chi2 = float(residuals.T @ precision @ residuals)
        rank = int(np.linalg.matrix_rank(design))
        return {
            "coefficients": coeffs,
            "covariance": covariance,
            "chi2": chi2,
            "rank": rank,
        }

    if sigma is None:
        raise ValueError("sigma is required for diagonal weighted fits")
    weights = 1.0 / sigma
    weighted_design = design * weights[:, None]
    weighted_y = y * weights
    coeffs, *_ = np.linalg.lstsq(weighted_design, weighted_y, rcond=None)
    covariance = np.linalg.pinv(weighted_design.T @ weighted_design)
    residuals = (y - design @ coeffs) / sigma
    chi2 = float(np.sum(residuals**2))
    rank = int(np.linalg.matrix_rank(weighted_design))
    return {
        "coefficients": coeffs,
        "covariance": covariance,
        "chi2": chi2,
        "rank": rank,
    }


def usable_bootstrap_matrix(bootstrap_y: np.ndarray | None, n_rows: int) -> np.ndarray | None:
    if bootstrap_y is None:
        return None
    matrix = np.asarray(bootstrap_y, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != n_rows or matrix.shape[1] < 2:
        return None
    return matrix


def fit_linear_model(
    rows: list[dict[str, Any]],
    design: np.ndarray,
    term_names: tuple[str, ...],
    bootstrap_y: np.ndarray | None = None,
) -> dict[str, Any]:
    if len(rows) < len(term_names):
        raise ValueError(
            f"The selected ansatz has {len(term_names)} coefficients; "
            f"need at least {len(term_names)} usable points, got {len(rows)}."
        )
    if design.shape != (len(rows), len(term_names)):
        raise ValueError(
            f"Internal fit-design shape mismatch: expected {(len(rows), len(term_names))}, got {design.shape}."
        )

    y = np.asarray([float(row["chi_hat"]) for row in rows], dtype=float)
    sigma = np.asarray([sigma_from_row(row) for row in rows], dtype=float)
    covariance_y: np.ndarray | None = None
    covariance_source = "diagonal_total_errors"
    bootstrap_section: dict[str, Any] | None = None
    bootstrap_matrix = usable_bootstrap_matrix(bootstrap_y, len(rows))

    if bootstrap_matrix is not None:
        stat_cov = block_statistical_covariance(rows, bootstrap_matrix)
        sys = np.asarray([finite_float(row.get("sys_err")) or 0.0 for row in rows], dtype=float)
        covariance_y = stat_cov + np.diag(sys**2)
        covariance_source = "block_bootstrap_statistical_covariance_plus_diagonal_interpolation_systematic"

    solved = solve_linear_ansatz(design, y, sigma=sigma, covariance_y=covariance_y)
    coeffs = np.asarray(solved["coefficients"], dtype=float)
    covariance = np.asarray(solved["covariance"], dtype=float)
    chi2 = float(solved["chi2"])
    rank = int(solved["rank"])
    if rank < len(term_names):
        raise ValueError(
            f"The selected ansatz is rank deficient for these points "
            f"(rank {rank} < {len(term_names)} coefficients)."
        )
    dof = int(max(0, len(rows) - rank))

    continuum_value = float(coeffs[0])
    continuum_err = float(math.sqrt(max(0.0, covariance[0, 0])))
    continuum_error_source = "linearized_parameter_covariance"
    coefficient_estimate_by_name = {
        name: float(value)
        for name, value in zip(term_names, coeffs, strict=True)
    }

    if bootstrap_matrix is not None:
        coefficient_samples: list[np.ndarray] = []
        for sample_index in range(bootstrap_matrix.shape[1]):
            y_sample = np.asarray(bootstrap_matrix[:, sample_index], dtype=float)
            if not np.all(np.isfinite(y_sample)):
                continue
            try:
                sample_solved = solve_linear_ansatz(design, y_sample, sigma=sigma, covariance_y=covariance_y)
            except Exception:
                continue
            coefficient_samples.append(np.asarray(sample_solved["coefficients"], dtype=float))
        if len(coefficient_samples) >= 2:
            coeff_sample_matrix = np.vstack(coefficient_samples)
            coeff_mean = np.nanmean(coeff_sample_matrix, axis=0)
            coeff_std = np.nanstd(coeff_sample_matrix, axis=0, ddof=1)
            coeff_p16 = np.nanpercentile(coeff_sample_matrix, 16, axis=0)
            coeff_p84 = np.nanpercentile(coeff_sample_matrix, 84, axis=0)
            continuum_value = float(coeff_mean[0])
            continuum_err = float(coeff_std[0])
            continuum_error_source = "bootstrap_refit_std"
            coefficient_estimate_by_name = {
                name: float(value)
                for name, value in zip(term_names, coeff_mean, strict=True)
            }
            bootstrap_section = {
                "samples_used": int(coeff_sample_matrix.shape[0]),
                "coefficient_mean_by_name": {
                    name: float(value)
                    for name, value in zip(term_names, coeff_mean, strict=True)
                },
                "coefficient_stat_err_by_name": {
                    name: float(value)
                    for name, value in zip(term_names, coeff_std, strict=True)
                },
                "coefficient_p16_by_name": {
                    name: float(value)
                    for name, value in zip(term_names, coeff_p16, strict=True)
                },
                "coefficient_p84_by_name": {
                    name: float(value)
                    for name, value in zip(term_names, coeff_p84, strict=True)
                },
                "continuum_mean": float(coeff_mean[0]),
                "continuum_stat_err": float(coeff_std[0]),
                "continuum_p16": float(coeff_p16[0]),
                "continuum_p84": float(coeff_p84[0]),
                "coefficient_samples": coeff_sample_matrix,
            }

    coeff_dict = {name: float(value) for name, value in zip(term_names, coeffs, strict=True)}
    err_dict = {
        name: float(math.sqrt(max(0.0, covariance[index, index])))
        for index, name in enumerate(term_names)
    }
    return {
        "term_names": list(term_names),
        "coefficients": [float(value) for value in coeffs],
        "coefficient_by_name": coeff_dict,
        "coefficient_estimate_by_name": coefficient_estimate_by_name,
        "coefficient_err_by_name": err_dict,
        "covariance": covariance,
        "data_covariance_source": covariance_source,
        "data_covariance": covariance_y,
        "rank": rank,
        "continuum_central_value": float(coeffs[0]),
        "continuum_central_err": float(math.sqrt(max(0.0, covariance[0, 0]))),
        "continuum_value": continuum_value,
        "continuum_err": continuum_err,
        "continuum_error_source": continuum_error_source,
        "bootstrap": bootstrap_section,
        "chi2": chi2,
        "dof": dof,
        "chi2_dof": float(chi2 / dof) if dof > 0 else None,
    }


def fit_combined_ansatz(rows: list[dict[str, Any]], bootstrap_y: np.ndarray | None = None) -> dict[str, Any]:
    if len(rows) < len(FIT_TERM_NAMES):
        raise ValueError(
            f"The combined Risch ansatz has {len(FIT_TERM_NAMES)} coefficients; "
            f"need at least {len(FIT_TERM_NAMES)} usable points, got {len(rows)}."
        )
    fit = fit_linear_model(rows, design_for_rows(rows), FIT_TERM_NAMES, bootstrap_y)
    fit.update(
        {
            "mode": "combined-risch",
            "formula": "chi_hat = c00 + c20*ahat2 + c40*ahat2^2 + c01*tau + c21*ahat2*tau + c02*tau^2",
            "fixed_smearing_formula": "tau=8*t_lat*ahat2",
            "n_points": int(len(rows)),
        }
    )
    return fit


def fit_separate_smearing_ansatz(
    rows: list[dict[str, Any]],
    *,
    order: int,
    bootstrap_y: np.ndarray | None = None,
) -> dict[str, Any]:
    term_names = separate_term_names(order)
    bootstrap_matrix = usable_bootstrap_matrix(bootstrap_y, len(rows))
    t_lat_values = sorted_unique([float(row_t_lat(row)) for row in rows if row_t_lat(row) is not None])
    smearing_fits: list[dict[str, Any]] = []
    failures: list[str] = []

    for t_lat in t_lat_values:
        row_indices = [
            index
            for index, row in enumerate(rows)
            if row_t_lat(row) is not None and abs(float(row_t_lat(row)) - t_lat) < 1e-10
        ]
        subset = [rows[index] for index in row_indices]
        if len(subset) < len(term_names):
            failures.append(
                f"{t_lat_label(t_lat)}: need {len(term_names)} usable points, got {len(subset)}"
            )
            continue
        ahat2 = np.asarray([float(row["ahat2"]) for row in subset], dtype=float)
        subset_bootstrap = bootstrap_matrix[row_indices, :] if bootstrap_matrix is not None else None
        try:
            subfit = fit_linear_model(
                subset,
                separate_design(ahat2, order),
                term_names,
                subset_bootstrap,
            )
        except ValueError as exc:
            failures.append(f"{t_lat_label(t_lat)}: {exc}")
            continue
        subfit.update(
            {
                "mode": "separate-linear" if order == 1 else "separate-quadratic",
                "order": int(order),
                "t_lat": float(t_lat),
                "eight_t_over_a2": float(t_lat_to_eight_t_over_a2(t_lat)),
                "t_over_a2_values": sorted_unique([float(row_t_lat(row)) for row in subset if row_t_lat(row) is not None]),
                "n_points": int(len(subset)),
                "formula": (
                    "chi_hat_tlat = c0 + c2*ahat2"
                    if order == 1
                    else "chi_hat_tlat = c0 + c2*ahat2 + c4*ahat2^2"
                ),
            }
        )
        smearing_fits.append(subfit)

    if failures:
        raise ValueError("Could not fit all fixed-smearing curves: " + "; ".join(failures))
    if not smearing_fits:
        raise ValueError("No fixed-smearing continuum fits could be built.")

    mode = "separate-linear" if order == 1 else "separate-quadratic"
    return {
        "mode": mode,
        "order": int(order),
        "term_names": list(term_names),
        "formula": (
            "For each fixed t_lat: chi_hat_tlat = c0^(t_lat) + c2^(t_lat)*ahat2"
            if order == 1
            else "For each fixed t_lat: chi_hat_tlat = c0^(t_lat) + c2^(t_lat)*ahat2 + c4^(t_lat)*ahat2^2"
        ),
        "smearing_fits": smearing_fits,
        "continuum_estimates": continuum_estimate_rows_from_smearing_fits(smearing_fits),
    }


def fit_continuum_ansatz(
    rows: list[dict[str, Any]],
    *,
    mode: str,
    bootstrap_y: np.ndarray | None = None,
) -> dict[str, Any]:
    if mode == "combined-risch":
        fit = fit_combined_ansatz(rows, bootstrap_y)
        fit["continuum_estimates"] = continuum_estimate_rows(fit)
        return fit
    return fit_separate_smearing_ansatz(rows, order=fit_mode_order(mode), bootstrap_y=bootstrap_y)


def predict_with_band(fit: dict[str, Any], design: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    coeffs = np.asarray(fit["coefficients"], dtype=float)
    covariance = np.asarray(fit["covariance"], dtype=float)
    y = design @ coeffs
    variance = np.einsum("ij,jk,ik->i", design, covariance, design)
    err = np.sqrt(np.maximum(0.0, variance))
    return y, err


def coefficient_sample_matrix(fit: dict[str, Any], n_terms: int) -> np.ndarray | None:
    raw_samples = (fit.get("bootstrap") or {}).get("coefficient_samples")
    if raw_samples is None:
        return None
    samples = np.asarray(raw_samples, dtype=float)
    if samples.ndim != 2 or samples.shape[0] < 2 or samples.shape[1] != n_terms:
        return None
    finite = np.all(np.isfinite(samples), axis=1)
    samples = samples[finite, :]
    return samples if samples.shape[0] >= 2 else None


def bootstrap_prediction_band(coefficient_samples: np.ndarray, design: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    y_samples = np.asarray(coefficient_samples, dtype=float) @ np.asarray(design, dtype=float).T
    finite = np.all(np.isfinite(y_samples), axis=1)
    y_samples = y_samples[finite, :]
    if y_samples.shape[0] < 2:
        raise ValueError("At least two finite bootstrap prediction samples are required.")
    return (
        np.nanmean(y_samples, axis=0),
        np.nanpercentile(y_samples, 16, axis=0),
        np.nanpercentile(y_samples, 84, axis=0),
    )


def prediction_band_for_fit(fit: dict[str, Any], design: np.ndarray) -> PredictionBand:
    samples = coefficient_sample_matrix(fit, np.asarray(design).shape[1])
    if samples is not None:
        try:
            y, y_low, y_high = bootstrap_prediction_band(samples, design)
            return PredictionBand(y=y, y_low=y_low, y_high=y_high, source="bootstrap_refit_percentile_16_84")
        except ValueError:
            pass
    y, err = predict_with_band(fit, design)
    return PredictionBand(y=y, y_low=y - err, y_high=y + err, source="linearized_parameter_covariance")


def fit_for_t_lat(fit: dict[str, Any], t_lat: float) -> dict[str, Any]:
    if fit.get("mode") != "combined-risch":
        for subfit in fit.get("smearing_fits", []):
            subfit_t_lat = finite_float(subfit.get("t_lat"))
            if subfit_t_lat is None:
                subfit_t_lat = finite_float((subfit.get("t_over_a2_values") or [None])[0])
            if subfit_t_lat is None:
                subfit_t_lat = finite_float(subfit.get("eight_t_over_a2"))
                subfit_t_lat = None if subfit_t_lat is None else subfit_t_lat / 8.0
            if subfit_t_lat is not None and abs(subfit_t_lat - float(t_lat)) < 1e-10:
                return subfit
        raise KeyError(f"No fit is available for {t_lat_label(t_lat)}.")
    return fit


def fit_for_smearing(fit: dict[str, Any], smearing_strength: float) -> dict[str, Any]:
    return fit_for_t_lat(fit, float(smearing_strength) / 8.0)


def design_grid_for_t_lat(fit: dict[str, Any], ahat2_grid: np.ndarray, t_lat: float) -> np.ndarray:
    if fit.get("mode") == "combined-risch":
        return smearing_design_grid(ahat2_grid, t_lat_to_eight_t_over_a2(t_lat))
    return separate_design(ahat2_grid, int(fit["order"]))


def design_grid_for_smearing(fit: dict[str, Any], ahat2_grid: np.ndarray, smearing_strength: float) -> np.ndarray:
    return design_grid_for_t_lat(fit, ahat2_grid, float(smearing_strength) / 8.0)


def continuum_estimate_rows_from_smearing_fits(smearing_fits: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for subfit in smearing_fits:
        rows.append(
            {
                "fit_mode": subfit.get("mode"),
                "t_lat": finite_float(subfit.get("t_lat")),
                "eight_t_over_a2": float(subfit["eight_t_over_a2"]),
                "t_over_a2_values": [float(value) for value in subfit.get("t_over_a2_values", [])],
                "continuum_value": finite_float(subfit.get("continuum_value")),
                "continuum_err": finite_float(subfit.get("continuum_err")),
                "continuum_error_source": subfit.get("continuum_error_source"),
                "chi2": finite_float(subfit.get("chi2")),
                "dof": subfit.get("dof"),
                "chi2_dof": finite_float(subfit.get("chi2_dof")),
                "n_points": int(subfit.get("n_points", 0)),
                "bootstrap_samples_used": (subfit.get("bootstrap") or {}).get("samples_used"),
            }
        )
    return rows


def continuum_estimate_rows(fit: dict[str, Any]) -> list[dict[str, Any]]:
    if fit.get("mode") != "combined-risch":
        return continuum_estimate_rows_from_smearing_fits(list(fit.get("smearing_fits", [])))
    return [
        {
            "fit_mode": fit.get("mode"),
            "t_lat": None,
            "eight_t_over_a2": None,
            "t_over_a2_values": None,
            "continuum_value": finite_float(fit.get("continuum_value")),
            "continuum_err": finite_float(fit.get("continuum_err")),
            "continuum_error_source": fit.get("continuum_error_source"),
            "chi2": finite_float(fit.get("chi2")),
            "dof": fit.get("dof"),
            "chi2_dof": finite_float(fit.get("chi2_dof")),
            "n_points": int(fit.get("n_points", 0)),
            "bootstrap_samples_used": (fit.get("bootstrap") or {}).get("samples_used"),
        }
    ]


def smearing_design_grid(ahat2_grid: np.ndarray, smearing_strength: float) -> np.ndarray:
    tau = smearing_strength * ahat2_grid
    return combined_design(ahat2_grid, tau)


def physical_tau_design_grid(ahat2_grid: np.ndarray, tau: float) -> np.ndarray:
    return combined_design(ahat2_grid, np.full_like(ahat2_grid, float(tau)))


def write_table(path: Path, rows: list[dict[str, Any]]) -> None:
    columns = [
        "beta",
        "epsilon_fixed",
        "epsilon_inferred",
        "scale_over_a",
        "t0_over_a2",
        "ahat2",
        "r_hat",
        "t_over_a2",
        "t_lat",
        "eight_t_over_a2",
        "tau",
        "chi_hat",
        "total_err",
        "stat_err",
        "sys_err",
        "used_in_fit",
        "fit_exclusion_reason",
        "fit_json",
        "analysis_dir",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# " + " ".join(columns) + "\n")
        for row in rows:
            out: list[str] = []
            for column in columns:
                value = row.get(column)
                if value is None:
                    out.append("nan")
                elif isinstance(value, float):
                    out.append(f"{value:.12g}")
                else:
                    out.append(str(value))
            handle.write(" ".join(out) + "\n")


def write_continuum_comparison_table(path: Path, rows: list[dict[str, Any]]) -> None:
    columns = [
        "fit_mode",
        "t_lat",
        "eight_t_over_a2",
        "continuum_value",
        "continuum_err",
        "continuum_error_source",
        "chi2",
        "dof",
        "chi2_dof",
        "n_points",
        "bootstrap_samples_used",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# " + " ".join(columns) + "\n")
        for row in rows:
            out: list[str] = []
            for column in columns:
                value = row.get(column)
                if value is None:
                    out.append("nan")
                elif isinstance(value, float):
                    out.append(f"{value:.12g}")
                else:
                    out.append(str(value))
            handle.write(" ".join(out) + "\n")


def sorted_unique(values: list[float], *, decimals: int = 12) -> list[float]:
    return sorted({round(float(value), decimals) for value in values})


def nonzero_flow_plot_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [row for row in rows if abs(float(row["t_over_a2"])) > 1e-12]


def color_for_index(index: int) -> str:
    palette = plt.get_cmap("tab10")
    r, g, b, _ = palette(index % 10)
    return f"#{int(255 * r):02x}{int(255 * g):02x}{int(255 * b):02x}"


def hex_to_rgba(color: str, alpha: float) -> str:
    color = color.lstrip("#")
    if len(color) != 6:
        return f"rgba(0,0,0,{alpha:g})"
    r = int(color[0:2], 16)
    g = int(color[2:4], 16)
    b = int(color[4:6], 16)
    return f"rgba({r},{g},{b},{alpha:g})"


def beta_marker_map(rows: list[dict[str, Any]]) -> dict[float | None, str]:
    markers = ["o", "s", "^", "D", "v", "P", "X", "<", ">"]
    betas = sorted({row.get("beta") for row in rows}, key=beta_sort_key)
    return {beta: markers[index % len(markers)] for index, beta in enumerate(betas)}


def save_continuum_png(
    path: Path,
    rows: list[dict[str, Any]],
    fit: dict[str, Any],
    *,
    scale: str,
    r_hat: float,
    eps_label: str,
    physical_tau: list[float],
) -> dict[str, Path]:
    fig, ax = plt.subplots(figsize=(6.2, 4.1), dpi=160)
    ahat2_values = [float(row["ahat2"]) for row in rows]
    x_max = max(ahat2_values) if ahat2_values else 0.0
    x_grid = np.linspace(0.0, max(x_max * 1.05, 1e-12), 240)
    t_lat_values = sorted_unique([float(row_t_lat(row)) for row in rows if row_t_lat(row) is not None])
    markers = beta_marker_map(rows)
    excluded_label_drawn = False

    for index, t_lat in enumerate(t_lat_values):
        color = color_for_index(index)
        subset = [row for row in rows if row_t_lat(row) is not None and abs(float(row_t_lat(row)) - t_lat) < 1e-10]
        subset_betas = sorted({row.get("beta") for row in subset}, key=beta_sort_key)
        used_betas = sorted({row.get("beta") for row in subset if bool(row.get("used_in_fit", True))}, key=beta_sort_key)
        first_used_beta = used_betas[0] if used_betas else (subset_betas[0] if subset_betas else None)
        for beta in subset_betas:
            beta_rows = [row for row in subset if row.get("beta") == beta]
            for used_in_fit, point_rows in (
                (True, [row for row in beta_rows if bool(row.get("used_in_fit", True))]),
                (False, [row for row in beta_rows if not bool(row.get("used_in_fit", True))]),
            ):
                if not point_rows:
                    continue
                x = np.asarray([float(row["ahat2"]) for row in point_rows], dtype=float)
                y = np.asarray([float(row["chi_hat"]) for row in point_rows], dtype=float)
                yerr = np.asarray([sigma_from_row(row) for row in point_rows], dtype=float)
                if used_in_fit:
                    ax.errorbar(
                        x,
                        y,
                        yerr=yerr,
                        fmt=markers[beta],
                        color=color,
                        ecolor=color,
                        capsize=2.0,
                        markersize=4.2,
                        linestyle="none",
                        label=t_lat_label(t_lat) if beta == first_used_beta else None,
                    )
                    label_color = color
                else:
                    ax.errorbar(
                        x,
                        y,
                        yerr=yerr,
                        fmt=markers[beta],
                        color="#7f7f7f",
                        ecolor="#7f7f7f",
                        markerfacecolor="none",
                        capsize=2.0,
                        markersize=4.2,
                        linestyle="none",
                        label="excluded from fit" if not excluded_label_drawn else None,
                    )
                    excluded_label_drawn = True
                    label_color = "#7f7f7f"
                for x_val, y_val, row in zip(x, y, point_rows, strict=True):
                    label = "" if row.get("beta") is None else f"{float(row['beta']):g}"
                    ax.annotate(label, (x_val, y_val), xytext=(3, 3), textcoords="offset points", fontsize=6, color=label_color)

        subfit = fit_for_t_lat(fit, t_lat)
        design = design_grid_for_t_lat(fit, x_grid, t_lat)
        band = prediction_band_for_fit(subfit, design)
        ax.plot(x_grid, band.y, color=color, linewidth=1.25)
        ax.fill_between(x_grid, band.y_low, band.y_high, color=color, alpha=0.12, linewidth=0)

        if fit.get("mode") != "combined-risch":
            ax.errorbar(
                [0.0],
                [float(subfit["continuum_value"])],
                yerr=[float(subfit["continuum_err"])],
                fmt="s",
                color=color,
                ecolor=color,
                capsize=2.0,
                markersize=5.0,
                label=rf"$a=0$, $t_{{\mathrm{{lat}}}}={t_lat:g}$",
            )

    physical_tau_values = physical_tau if fit.get("mode") == "combined-risch" else []
    for tau in physical_tau_values:
        design = physical_tau_design_grid(x_grid, tau)
        band = prediction_band_for_fit(fit, design)
        ax.plot(x_grid, band.y, color="#444444", linewidth=1.0, linestyle="--", label=rf"$\tau={tau:g}$")
        ax.fill_between(x_grid, band.y_low, band.y_high, color="#444444", alpha=0.07, linewidth=0)

    if fit.get("mode") == "combined-risch":
        ax.errorbar(
            [0.0],
            [float(fit["continuum_value"])],
            yerr=[float(fit["continuum_err"])],
            fmt="s",
            color="#d62728",
            ecolor="#d62728",
            capsize=2.0,
            markersize=5.0,
            label=r"$a=0$",
        )
    ax.set_xlabel(r"$(a/r_0)^2$" if scale == "r0" else r"$a^2/(8t_0)$")
    ax.set_ylabel(r"$\hat{\chi}(\hat{r})$")
    ax.set_title(f"r_hat={r_hat:g}, {eps_label}")
    ax.set_xlim(left=0.0)
    ax.tick_params(direction="in", top=True, right=True)
    ax.legend(frameon=True, fontsize=7)
    fig.tight_layout()
    written = save_matplotlib_figure(fig, path)
    plt.close(fig)
    if not written:
        raise RuntimeError(f"No continuum plot exports could be written for {path}.")
    return written


def add_plotly_band(
    fig: Any,
    x: np.ndarray,
    y: np.ndarray,
    y_low: np.ndarray,
    y_high: np.ndarray,
    *,
    color: str,
    name: str,
    dash: str | None = None,
) -> None:
    go = _get_plotly()
    fig.add_trace(
        go.Scatter(
            x=np.concatenate([x, x[::-1]]),
            y=np.concatenate([y_high, y_low[::-1]]),
            fill="toself",
            fillcolor=hex_to_rgba(color, 0.16),
            line={"color": "rgba(0,0,0,0)"},
            hoverinfo="skip",
            showlegend=False,
            name=f"{name} band",
        )
    )
    line: dict[str, Any] = {"color": color, "width": 2}
    if dash is not None:
        line["dash"] = dash
    fig.add_trace(go.Scatter(x=x, y=y, mode="lines", line=line, name=name))


def save_continuum_html(
    path: Path,
    rows: list[dict[str, Any]],
    fit: dict[str, Any],
    *,
    scale: str,
    r_hat: float,
    eps_label: str,
    physical_tau: list[float],
) -> None:
    go = _get_plotly()
    x_hover_label = "(a/r<sub>0</sub>)<sup>2</sup>" if scale == "r0" else "a<sup>2</sup>/(8 t<sub>0</sub>)"
    fig = go.Figure()
    ahat2_values = [float(row["ahat2"]) for row in rows]
    x_max = max(ahat2_values) if ahat2_values else 0.0
    x_grid = np.linspace(0.0, max(x_max * 1.05, 1e-12), 240)
    t_lat_values = sorted_unique([float(row_t_lat(row)) for row in rows if row_t_lat(row) is not None])

    for index, t_lat in enumerate(t_lat_values):
        color = color_for_index(index)
        subset = [row for row in rows if row_t_lat(row) is not None and abs(float(row_t_lat(row)) - t_lat) < 1e-10]
        for used_in_fit, point_rows in (
            (True, [row for row in subset if bool(row.get("used_in_fit", True))]),
            (False, [row for row in subset if not bool(row.get("used_in_fit", True))]),
        ):
            if not point_rows:
                continue
            marker = {"color": color, "size": 7}
            name = t_lat_legend_label(t_lat, "data")
            if not used_in_fit:
                marker = {"color": "#7f7f7f", "size": 7, "symbol": "circle-open"}
                name = t_lat_legend_label(t_lat, "excluded")
            fig.add_trace(
                go.Scatter(
                    x=[float(row["ahat2"]) for row in point_rows],
                    y=[float(row["chi_hat"]) for row in point_rows],
                    error_y={"type": "data", "array": [sigma_from_row(row) for row in point_rows], "visible": True},
                    mode="markers+text",
                    text=[f"{float(row['beta']):g}" if row.get("beta") is not None else "" for row in point_rows],
                    textposition="top center",
                    marker=marker,
                    name=name,
                    customdata=[
                        [
                            row.get("beta"),
                            row.get("epsilon_fixed"),
                            row.get("epsilon_inferred"),
                            row.get("t0_over_a2"),
                            row.get("t_over_a2"),
                            row.get("eight_t_over_a2"),
                            row.get("tau"),
                            row.get("stat_err"),
                            row.get("sys_err"),
                            row.get("total_err"),
                            row.get("analysis_dir"),
                            row.get("used_in_fit", True),
                            row.get("fit_exclusion_reason"),
                        ]
                        for row in point_rows
                    ],
                    hovertemplate=(
                        "beta=%{customdata[0]}<br>"
                        "fixed &epsilon;<sub>1</sub>=%{customdata[1]}<br>"
                        "inferred &epsilon;<sub>1</sub>=%{customdata[2]}<br>"
                        "t0/a^2=%{customdata[3]:g}<br>"
                        "t_lat=%{customdata[4]:g}<br>"
                        "tau=%{customdata[6]:g}<br>"
                        f"{x_hover_label}=%{{x:g}}<br>"
                        "&chi;&#770;=%{y:g}<br>"
                        "stat=%{customdata[7]:g}<br>"
                        "sys=%{customdata[8]:g}<br>"
                        "total=%{customdata[9]:g}<br>"
                        "used in fit=%{customdata[11]}<br>"
                        "fit exclusion=%{customdata[12]}<br>"
                        "%{customdata[10]}<extra></extra>"
                    ),
                )
            )

        subfit = fit_for_t_lat(fit, t_lat)
        design = design_grid_for_t_lat(fit, x_grid, t_lat)
        band = prediction_band_for_fit(subfit, design)
        add_plotly_band(fig, x_grid, band.y, band.y_low, band.y_high, color=color, name=t_lat_legend_label(t_lat, "fit"))

        if fit.get("mode") != "combined-risch":
            fig.add_trace(
                go.Scatter(
                    x=[0.0],
                    y=[float(subfit["continuum_value"])],
                    error_y={"type": "data", "array": [float(subfit["continuum_err"])], "visible": True},
                    mode="markers",
                    marker={"symbol": "square", "color": color, "size": 9},
                    name=continuum_t_lat_legend_label(t_lat),
                    hovertemplate=(
                        f"{x_hover_label}=0<br>"
                        "t_lat=%{customdata[0]:g}<br>"
                        "&chi;&#770;=%{y:g}<br>"
                        "err=%{customdata[1]:g}<br>"
                        "chi2/dof=%{customdata[2]}<extra>a=0</extra>"
                    ),
                    customdata=[
                        [
                            float(t_lat),
                            float(subfit["continuum_err"]),
                            subfit.get("chi2_dof"),
                        ]
                    ],
                )
            )

    physical_tau_values = physical_tau if fit.get("mode") == "combined-risch" else []
    for tau in physical_tau_values:
        design = physical_tau_design_grid(x_grid, tau)
        band = prediction_band_for_fit(fit, design)
        add_plotly_band(fig, x_grid, band.y, band.y_low, band.y_high, color="#444444", name=f"physical tau={tau:g}", dash="dash")

    if fit.get("mode") == "combined-risch":
        fig.add_trace(
            go.Scatter(
                x=[0.0],
                y=[float(fit["continuum_value"])],
                error_y={"type": "data", "array": [float(fit["continuum_err"])], "visible": True},
                mode="markers",
                marker={"symbol": "square", "color": "#d62728", "size": 9},
                name="a=0, shared continuum",
                hovertemplate=f"{x_hover_label}=0<br>&chi;&#770;=%{{y:g}}<extra>a=0</extra>",
            )
        )
    fig.update_layout(
        title=f"Continuum fit at r_hat={r_hat:g}, {eps_label}",
        template="plotly_white",
        xaxis_title="$(a/r_0)^2$" if scale == "r0" else "$a^2/(8t_0)$",
        yaxis_title="$\\hat{\\chi}(\\hat{r})$",
        margin={"l": 58, "r": 12, "t": 48, "b": 82},
    )
    fig.update_xaxes(rangemode="tozero", automargin=True, title_standoff=16)
    _write_figure_html(fig, path)


def relative_uncertainty_entries(
    inputs: list[ContinuumInput],
    flow_times: list[float],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    entries: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    for item in inputs:
        for flow_time in flow_times:
            path = item.flow_paths.get(flow_time)
            if path is None:
                skipped.append({"analysis_dir": str(item.analysis_dir), "t_over_a2": flow_time, "reason": "missing flow fit"})
                continue
            fit = load_json(path)
            grid = fit.get("grid") or {}
            x = np.asarray(grid.get("x", []), dtype=float)
            average = np.asarray(grid.get("average", []), dtype=float)
            total_rel_raw = grid.get("total_rel")
            if total_rel_raw is not None:
                total_rel = np.asarray(total_rel_raw, dtype=float)
            else:
                total_err = np.asarray(grid.get("total_err", []), dtype=float)
                denominator = np.where(np.abs(average) > 0, np.abs(average), np.nan)
                total_rel = total_err / denominator
            if x.size == 0 or total_rel.size != x.size:
                skipped.append({"analysis_dir": str(item.analysis_dir), "t_over_a2": flow_time, "reason": "missing relative uncertainty grid"})
                continue
            valid = np.isfinite(x) & np.isfinite(total_rel)
            if np.count_nonzero(valid) < 2:
                skipped.append({"analysis_dir": str(item.analysis_dir), "t_over_a2": flow_time, "reason": "non-finite relative uncertainty grid"})
                continue
            order = np.argsort(x[valid])
            entries.append(
                {
                    "analysis_dir": str(item.analysis_dir),
                    "fit_json": str(path),
                    "beta": item.beta,
                    "epsilon_inferred": item.inferred_epsilon,
                    "t_over_a2": float(flow_time),
                    "t_lat": float(flow_time),
                    "eight_t_over_a2": float(8.0 * flow_time),
                    "x": x[valid][order],
                    "total_rel": total_rel[valid][order],
                    "x_min": float(np.min(x[valid])),
                    "x_max": float(np.max(x[valid])),
                }
            )
    return entries, skipped


def common_rhat_range(entries: list[dict[str, Any]]) -> tuple[float | None, float | None]:
    if not entries:
        return None, None
    lo = max(float(entry["x_min"]) for entry in entries)
    hi = min(float(entry["x_max"]) for entry in entries)
    if lo >= hi:
        return None, None
    return lo, hi


def use_log_relative_axis(entries: list[dict[str, Any]]) -> bool:
    values = np.concatenate([np.asarray(entry["total_rel"], dtype=float) for entry in entries]) if entries else np.asarray([])
    positive = values[np.isfinite(values) & (values > 0)]
    if positive.size < 2:
        return False
    return float(np.nanmax(positive) / np.nanmin(positive)) > 100.0


def save_rhat_guide_png(
    path: Path,
    entries: list[dict[str, Any]],
    *,
    eps_label: str,
    scale: str,
    common_range: tuple[float | None, float | None],
) -> dict[str, Path]:
    betas = sorted({entry.get("beta") for entry in entries}, key=beta_sort_key)
    if not betas:
        raise ValueError("No relative-uncertainty entries are available.")
    fig, axes = plt.subplots(len(betas), 1, figsize=(6.4, max(2.2, 1.9 * len(betas))), dpi=160, sharex=True)
    if len(betas) == 1:
        axes = [axes]

    flow_values = sorted_unique([float(entry["t_lat"]) for entry in entries])
    flow_color = {flow: color_for_index(index) for index, flow in enumerate(flow_values)}
    x_lo, x_hi = common_range
    for axis, beta in zip(axes, betas, strict=True):
        beta_entries = [entry for entry in entries if entry.get("beta") == beta]
        for entry in sorted(beta_entries, key=lambda row: float(row["t_lat"])):
            flow = float(entry["t_lat"])
            axis.plot(
                np.asarray(entry["x"], dtype=float),
                np.asarray(entry["total_rel"], dtype=float),
                color=flow_color[flow],
                linewidth=1.1,
                label=t_lat_label(flow),
            )
        if x_lo is not None and x_hi is not None:
            axis.axvspan(float(x_lo), float(x_hi), color="0.5", alpha=0.08, linewidth=0)
            axis.axvline(float(x_lo), color="0.55", linewidth=0.7, linestyle=":")
            axis.axvline(float(x_hi), color="0.55", linewidth=0.7, linestyle=":")
        if use_log_relative_axis(beta_entries):
            positive = np.concatenate([np.asarray(entry["total_rel"], dtype=float) for entry in beta_entries])
            positive = positive[np.isfinite(positive) & (positive > 0)]
            axis.set_yscale("log")
            axis.set_ylim(bottom=max(float(np.nanmin(positive)) * 0.7, 1e-8))
        else:
            axis.set_ylim(bottom=0.0)
        title = "beta=unknown" if beta is None else f"beta={beta:g}"
        axis.set_title(title, fontsize=9)
        axis.set_ylabel(r"$\Delta_{\rm tot}/\hat\chi$")
        axis.tick_params(direction="in", top=True, right=True)
        axis.legend(frameon=True, fontsize=7, ncols=min(3, max(1, len(flow_values))))

    axes[-1].set_xlabel("r / r0" if scale == "r0" else "r / sqrt(8 t0)")
    fig.suptitle(f"Relative interpolation uncertainty, {eps_label}", fontsize=10)
    fig.tight_layout()
    written = save_matplotlib_figure(fig, path)
    plt.close(fig)
    if not written:
        raise RuntimeError(f"No r-hat guide plot exports could be written for {path}.")
    return written


def save_rhat_guide_html(
    path: Path,
    entries: list[dict[str, Any]],
    *,
    eps_label: str,
    scale: str,
    common_range: tuple[float | None, float | None],
) -> None:
    go = _get_plotly()
    from plotly.subplots import make_subplots

    betas = sorted({entry.get("beta") for entry in entries}, key=beta_sort_key)
    if not betas:
        raise ValueError("No relative-uncertainty entries are available.")
    fig = make_subplots(
        rows=len(betas),
        cols=1,
        shared_xaxes=True,
        subplot_titles=["beta=unknown" if beta is None else f"beta={beta:g}" for beta in betas],
    )
    flow_values = sorted_unique([float(entry["t_lat"]) for entry in entries])
    flow_color = {flow: color_for_index(index) for index, flow in enumerate(flow_values)}
    x_lo, x_hi = common_range
    for row_index, beta in enumerate(betas, start=1):
        beta_entries = [entry for entry in entries if entry.get("beta") == beta]
        for entry in sorted(beta_entries, key=lambda row: float(row["t_lat"])):
            flow = float(entry["t_lat"])
            fig.add_trace(
                go.Scatter(
                    x=[float(value) for value in entry["x"]],
                    y=[float(value) for value in entry["total_rel"]],
                    mode="lines",
                    line={"color": flow_color[flow], "width": 1.7},
                    name=t_lat_label(flow),
                    legendgroup=f"flow-{flow:g}",
                    showlegend=row_index == 1,
                    customdata=[
                        [entry.get("beta"), entry.get("epsilon_inferred"), entry.get("t_over_a2"), entry.get("fit_json")]
                        for _ in entry["x"]
                    ],
                    hovertemplate=(
                        "beta=%{customdata[0]}<br>"
                        "&epsilon;<sub>1</sub>=%{customdata[1]}<br>"
                        "t_lat=%{customdata[2]:g}<br>"
                        "r_hat=%{x:g}<br>"
                        "Delta_tot/chi=%{y:g}<br>"
                        "%{customdata[3]}<extra></extra>"
                    ),
                ),
                row=row_index,
                col=1,
            )
        if x_lo is not None and x_hi is not None:
            fig.add_vrect(
                x0=float(x_lo),
                x1=float(x_hi),
                fillcolor="gray",
                opacity=0.08,
                line_width=0,
                row=row_index,
                col=1,
            )
        if use_log_relative_axis(beta_entries):
            fig.update_yaxes(type="log", row=row_index, col=1)

    fig.update_layout(
        title=f"Relative interpolation uncertainty, {eps_label}",
        template="plotly_white",
        yaxis_title="Delta_tot / chi_hat",
    )
    fig.update_xaxes(title_text="r / r0" if scale == "r0" else "r / sqrt(8 t0)", row=len(betas), col=1)
    for row_index, beta in enumerate(betas, start=1):
        beta_entries = [entry for entry in entries if entry.get("beta") == beta]
        axis_options: dict[str, Any] = {"title_text": "Delta_tot / chi_hat"}
        if not use_log_relative_axis(beta_entries):
            axis_options["rangemode"] = "tozero"
        fig.update_yaxes(row=row_index, col=1, **axis_options)
    _write_figure_html(fig, path)


def save_rhat_guide_outputs(
    output_dir: Path,
    entries: list[dict[str, Any]],
    skipped: list[dict[str, Any]],
    *,
    eps_label: str,
    scale: str,
    betas: list[float],
    flow_times: list[float],
) -> dict[str, Path]:
    common_range = common_rhat_range(entries)
    beta_token = "betas_" + "_".join(filename_token(f"{value:g}") for value in betas)
    flow_token = "s_" + "_".join(filename_token(f"{8.0 * value:g}") for value in flow_times)
    token = f"rhat_guide__{beta_token}__{flow_token}__{scale}"
    png_path = output_dir / f"{token}.png"
    html_path = output_dir / f"{token}.html"
    summary_path = output_dir / f"{token}.json"
    plot_paths = save_rhat_guide_png(png_path, entries, eps_label=eps_label, scale=scale, common_range=common_range)
    save_rhat_guide_html(html_path, entries, eps_label=eps_label, scale=scale, common_range=common_range)
    save_json(
        summary_path,
        {
            "schema_version": 1,
            "scale": scale,
            "betas": [float(value) for value in betas],
            "t_over_a2_values": [float(value) for value in flow_times],
            "t_lat_values": [float(value) for value in flow_times],
            "eight_t_over_a2_values": [float(8.0 * value) for value in flow_times],
            "eps_label": eps_label,
            "common_r_hat_range": {
                "min": common_range[0],
                "max": common_range[1],
            },
            "entries": [
                {
                    key: value
                    for key, value in entry.items()
                    if key not in {"x", "total_rel"}
                }
                for entry in entries
            ],
            "skipped": skipped,
            "paths": {
                "plot_png": str(plot_paths.get("png", png_path)),
                "plot_pdf": str(plot_paths["pdf"]) if "pdf" in plot_paths else None,
                "plot_html": str(html_path),
                "summary": str(summary_path),
                "plot_formats": {fmt: str(fmt_path) for fmt, fmt_path in plot_paths.items()},
            },
            "saved_at": datetime.now(timezone.utc).isoformat(),
        },
    )
    return {**plot_paths, "html": html_path, "summary": summary_path}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Continuum fit from precomputed Creutz-ratio interpolation fits."
    )
    parser.add_argument(
        "paths",
        nargs="*",
        type=Path,
        default=[DEFAULT_INPUT_ROOT],
        help=f"creutz_ratio_analysis dirs or roots containing them. Default: {DEFAULT_INPUT_ROOT}.",
    )
    parser.add_argument("-o", "--output-dir", type=Path, default=DEFAULT_OUTPUT_ROOT, help=f"Output directory. Default: {DEFAULT_OUTPUT_ROOT}.")
    parser.add_argument("--scale", choices=["r0", "t0"], default="r0", help="Dimensionless Creutz-ratio scale to read from fits. Default: r0.")
    parser.add_argument("--r-hat", type=float, help="Target r_hat where the average interpolation fit is evaluated.")
    parser.add_argument("--betas", nargs="+", type=float, default=list(DEFAULT_BETAS), help="Beta values to use. Default: 2.3 2.4 2.5 2.6.")
    parser.add_argument(
        "--exclude-fit-beta",
        "--exclude-beta-from-fit",
        dest="exclude_fit_beta",
        nargs="+",
        type=float,
        default=[],
        metavar="BETA",
        help="Beta values to keep in tables/plots but exclude from the continuum fit.",
    )
    parser.add_argument(
        "--eps1",
        nargs="+",
        type=float,
        default=list(DEFAULT_EPS1),
        help="epsilon1 values. Give one value for all betas, or one per beta in --betas order. Default: 0.",
    )
    parser.add_argument(
        "--flow-time",
        "--t-lat",
        dest="flow_time",
        nargs="+",
        type=float,
        help="One or more Wilson-loop t_lat values to include. Default: 0 0.25 0.5 0.75 1.0.",
    )
    parser.add_argument(
        "--eight-t-over-a2",
        nargs="+",
        type=float,
        help="Legacy alias accepting 8*t_lat values; prefer --t-lat.",
    )
    parser.add_argument(
        "--hide-zero-flow-in-plot",
        action="store_true",
        help="Do not draw t_lat=0 in the continuum PNG/HTML plots; fits and tables still include it.",
    )
    parser.add_argument("--physical-tau", nargs="*", type=float, default=[], help="Optional dashed physical-flow curves at fixed tau=t_lat/(t0/a^2).")
    parser.add_argument(
        "--fit-mode",
        choices=list(FIT_MODES),
        default="separate-linear",
        help=(
            "Continuum ansatz. separate-linear/quadratic fits each fixed t_lat independently; "
            "combined-risch uses the shared six-parameter Risch ansatz. Default: separate-linear."
        ),
    )
    parser.add_argument("--eps-map", nargs="+", metavar="BETA=EPS", help="Per-beta fixed-epsilon values used to select the trajectory. Overrides --eps1.")
    parser.add_argument("--eps-label", help="Label for the fixed-epsilon trajectory in plots.")
    parser.add_argument("--eps-tol", type=float, default=1e-8, help="Tolerance for matching inferred epsilon1 to --eps-map. Default: 1e-8.")
    parser.add_argument("--no-eps-prompt", action="store_true", help="Do not prompt for ambiguous fixed-epsilon choices.")
    parser.add_argument("--t0-map", nargs="+", metavar="BETA=T0_OVER_A2", help="Per-beta t0/a^2 values for ahat^2 and tau. Overrides built-in defaults for beta 2.3-2.6.")
    parser.add_argument("--no-t0-prompt", action="store_true", help="Do not prompt for missing t0/a^2 values.")
    parser.add_argument("--continuum-bootstrap-samples", type=int, default=200, help="Maximum bootstrap replicas used in the continuum fit. Default: 200.")
    parser.add_argument("--continuum-bootstrap-seed", type=int, default=24680, help="Seed for Gaussian fallback replicas. Default: 24680.")
    parser.add_argument("--no-continuum-bootstrap", action="store_true", help="Use only diagonal total errors in the continuum fit.")
    parser.add_argument("--list", action="store_true", help="List discovered inputs and common t_lat values, then exit.")
    return parser


def validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    if args.r_hat is not None and args.r_hat <= 0:
        parser.error("--r-hat must be positive")
    if not args.betas:
        parser.error("--betas must contain at least one value")
    for beta in args.betas:
        if not np.isfinite(float(beta)):
            parser.error("--betas values must be finite")
    selected_betas = sorted_unique([float(beta) for beta in args.betas])
    excluded_betas = sorted_unique([float(beta) for beta in args.exclude_fit_beta])
    for beta in excluded_betas:
        if not np.isfinite(float(beta)):
            parser.error("--exclude-fit-beta values must be finite")
        if not float_in_list(beta, args.betas):
            parser.error(f"--exclude-fit-beta {beta:g} is not present in --betas")
    if len(excluded_betas) >= len(selected_betas):
        parser.error("--exclude-fit-beta cannot exclude every selected beta")
    if not args.eps1 and not args.eps_map:
        parser.error("--eps1 must contain at least one value unless --eps-map is given")
    if args.flow_time is not None and args.eight_t_over_a2 is not None:
        parser.error("Use only one of --flow-time and --eight-t-over-a2")
    if args.eps_tol <= 0:
        parser.error("--eps-tol must be positive")
    if args.continuum_bootstrap_samples < 2:
        parser.error("--continuum-bootstrap-samples must be >= 2")
    for tau in args.physical_tau:
        if tau < 0:
            parser.error("--physical-tau values must be non-negative")
    if args.physical_tau and args.fit_mode != "combined-risch":
        parser.error("--physical-tau curves are only available with --fit-mode combined-risch")


def make_eps_label(args: argparse.Namespace, eps_info: dict[str, Any]) -> str:
    if args.eps_label:
        return args.eps_label
    values = sorted({value for value in (eps_info.get("by_beta") or {}).values() if value is not None})
    if not values:
        return r"fixed $\epsilon_1$"
    return r"fixed $\epsilon_1=$" + ",".join(f"{value:g}" for value in values)


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    validate_args(parser, args)

    try:
        betas = [float(value) for value in args.betas]
        analysis_dirs = discover_analysis_dirs(args.paths)
        inputs = [item for path in analysis_dirs if (item := load_continuum_input(path, args.scale)) is not None]
        inputs = filter_inputs_by_beta(inputs, betas)
        eps_map = resolve_eps_map_from_args(args, betas)
        t0_map = resolve_t0_map_from_args(args, betas)
    except (FileNotFoundError, ValueError) as exc:
        parser.error(str(exc))

    if not inputs:
        beta_text = ", ".join(f"{value:g}" for value in args.betas)
        parser.error(f"No creutz_ratio_analysis inputs with fits_{args.scale} were found for beta values: {beta_text}.")

    if args.list:
        try:
            list_inputs, _list_eps_by_dir, eps_info, eps_skipped = resolve_fixed_eps_selection(
                inputs,
                eps_map,
                prompt=False,
                eps_tol=args.eps_tol,
            )
        except ValueError as exc:
            parser.error(str(exc))
        common = common_flow_times(list_inputs)
        print("Common t_lat: " + (", ".join(f"{value:g}" for value in common) if common else "none"))
        print("Selected fixed eps: " + json.dumps(eps_info["by_beta"], sort_keys=True))
        for item in sorted(list_inputs, key=lambda row: (beta_sort_key(row.beta), str(row.analysis_dir))):
            t0_text = "none" if item.summary_t0_over_a2 is None else f"{item.summary_t0_over_a2:g}"
            scale_text = "none" if item.scale_over_a is None else f"{item.scale_over_a:g}"
            print(
                f"beta={item.beta}, epsilon1={item.inferred_epsilon}, "
                f"scale_over_a={scale_text}, summary_t0/a^2={t0_text}: {item.analysis_dir}"
            )
        if eps_skipped:
            print(f"Skipped {len(eps_skipped)} input(s) due to epsilon selection.")
        return 0

    try:
        selected_inputs, fixed_eps_by_dir, eps_info, eps_skipped = resolve_fixed_eps_selection(
            inputs,
            eps_map,
            prompt=not args.no_eps_prompt,
            eps_tol=args.eps_tol,
        )
        if not selected_inputs:
            raise ValueError("No inputs remain after fixed-epsilon selection.")
        flow_times = resolve_flow_times(args, selected_inputs)
    except ValueError as exc:
        parser.error(str(exc))

    eps_label = make_eps_label(args, eps_info)
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.r_hat is None:
        guide_entries, guide_skipped = relative_uncertainty_entries(selected_inputs, flow_times)
        guide_skipped.extend(eps_skipped)
        if not guide_entries:
            parser.error("No relative-uncertainty grids were available for the selected beta/eps/flow inputs.")
        try:
            guide_paths = save_rhat_guide_outputs(
                output_dir,
                guide_entries,
                guide_skipped,
                eps_label=eps_label,
                scale=args.scale,
                betas=betas,
                flow_times=flow_times,
            )
        except ValueError as exc:
            parser.error(str(exc))
        common_lo, common_hi = common_rhat_range(guide_entries)
        print(f"Saved r_hat guide plot: {guide_paths['html']}")
        print(f"Saved r_hat guide summary: {guide_paths['summary']}")
        if common_lo is not None and common_hi is not None:
            print(f"Common r_hat range across selected beta/flow inputs: {common_lo:.6g} to {common_hi:.6g}")
        print("Run again with --r-hat VALUE to perform the continuum fit.")
        return 0

    try:
        t0_by_dir, t0_info = resolve_t0_values(selected_inputs, t0_map, prompt=not args.no_t0_prompt)
        r_hat = resolve_r_hat(args)
    except ValueError as exc:
        parser.error(str(exc))

    rows, skipped = build_rows(
        selected_inputs,
        scale=args.scale,
        r_hat=r_hat,
        flow_times=flow_times,
        fixed_eps_by_dir=fixed_eps_by_dir,
        t0_by_dir=t0_by_dir,
    )
    skipped.extend(eps_skipped)
    if not rows:
        parser.error("No usable continuum points remain after interpolation to the selected r_hat.")
    fit_excluded_betas = sorted_unique([float(value) for value in args.exclude_fit_beta])
    fit_rows = mark_fit_usage(rows, fit_excluded_betas)
    if not fit_rows:
        parser.error("No usable continuum points remain after applying --exclude-fit-beta.")
    if args.fit_mode == "combined-risch" and len(fit_rows) < len(FIT_TERM_NAMES):
        parser.error(f"Need at least {len(FIT_TERM_NAMES)} usable points for the combined Risch ansatz; got {len(fit_rows)}.")

    bootstrap_matrix: np.ndarray | None = None
    bootstrap_metadata: list[dict[str, Any]] = []
    if not args.no_continuum_bootstrap:
        inputs_by_dir = {str(item.analysis_dir): item for item in selected_inputs}
        bootstrap_matrix, bootstrap_metadata = build_interpolated_bootstrap_matrix(
            fit_rows,
            inputs_by_dir,
            scale=args.scale,
            max_samples=args.continuum_bootstrap_samples,
            seed=args.continuum_bootstrap_seed,
        )

    try:
        fit = fit_continuum_ansatz(fit_rows, mode=args.fit_mode, bootstrap_y=bootstrap_matrix)
    except ValueError as exc:
        parser.error(str(exc))

    s_token = "s_" + "_".join(filename_token(f"{8.0 * value:g}") for value in flow_times)
    mode_token = args.fit_mode.replace("-", "_")
    token = f"rhat_{filename_token(f'{r_hat:g}')}__{s_token}__{args.scale}__{mode_token}"
    if fit_excluded_betas:
        excluded_token = "_".join(filename_token(f"{value:g}") for value in fit_excluded_betas)
        token += f"__fit_excludes_beta_{excluded_token}"
    plot_rows = nonzero_flow_plot_rows(rows) if args.hide_zero_flow_in_plot else rows
    if not plot_rows:
        parser.error("--hide-zero-flow-in-plot removed all continuum points; include at least one nonzero t_lat value.")
    plot_token = f"{token}__plot_no_s0" if args.hide_zero_flow_in_plot else token

    table_path = output_dir / f"continuum_points__{token}.dat"
    comparison_path = output_dir / f"continuum_limits__{token}.dat"
    summary_path = output_dir / f"continuum_summary__{token}.json"
    png_path = output_dir / f"continuum_plot__{plot_token}.png"
    html_path = output_dir / f"continuum_plot__{plot_token}.html"
    bootstrap_path = output_dir / f"continuum_bootstrap__{token}.npz"

    continuum_estimates = continuum_estimate_rows(fit)
    write_table(table_path, rows)
    write_continuum_comparison_table(comparison_path, continuum_estimates)
    if bootstrap_matrix is not None:
        npz_payload: dict[str, Any] = {
            "y_bootstrap": bootstrap_matrix,
            "row_analysis_dir": np.asarray([str(row["analysis_dir"]) for row in fit_rows]),
            "row_beta": np.asarray([float(row["beta"]) if row.get("beta") is not None else np.nan for row in fit_rows], dtype=float),
            "row_t_over_a2": np.asarray([float(row["t_over_a2"]) for row in fit_rows], dtype=float),
            "row_t_lat": np.asarray(
                [float(value) if (value := row_t_lat(row)) is not None else np.nan for row in fit_rows],
                dtype=float,
            ),
            "row_ahat2": np.asarray([float(row["ahat2"]) for row in fit_rows], dtype=float),
            "row_tau": np.asarray([float(row["tau"]) for row in fit_rows], dtype=float),
        }
        if fit.get("mode") == "combined-risch":
            raw_coefficient_samples = (fit.get("bootstrap") or {}).get("coefficient_samples")
            npz_payload["coefficient_samples"] = np.asarray(
                raw_coefficient_samples if raw_coefficient_samples is not None else [],
                dtype=float,
            )
            npz_payload["term_names"] = np.asarray(FIT_TERM_NAMES)
        else:
            for subfit in fit.get("smearing_fits", []):
                suffix = filename_token(f"{float(subfit['eight_t_over_a2']):g}")
                raw_coefficient_samples = (subfit.get("bootstrap") or {}).get("coefficient_samples")
                npz_payload[f"coefficient_samples_s_{suffix}"] = np.asarray(
                    raw_coefficient_samples if raw_coefficient_samples is not None else [],
                    dtype=float,
                )
                npz_payload[f"term_names_s_{suffix}"] = np.asarray(subfit.get("term_names", []))
        np.savez(bootstrap_path, **npz_payload)
    continuum_plot_paths = save_continuum_png(
        png_path,
        plot_rows,
        fit,
        scale=args.scale,
        r_hat=r_hat,
        eps_label=eps_label,
        physical_tau=args.physical_tau,
    )
    save_continuum_html(
        html_path,
        plot_rows,
        fit,
        scale=args.scale,
        r_hat=r_hat,
        eps_label=eps_label,
        physical_tau=args.physical_tau,
    )
    x_axis_definition = (
        "ahat2 = (a/r0)^2 = 1/(r0_over_a)^2"
        if args.scale == "r0"
        else "ahat2 = a^2/(8 t0) = 1/(8 t0_over_a2)"
    )
    save_json(
        summary_path,
        {
            "schema_version": 3,
            "scale": args.scale,
            "r_hat": float(r_hat),
            "t_over_a2_values": [float(value) for value in flow_times],
            "t_lat_values": [float(value) for value in flow_times],
            "eight_t_over_a2_values": [float(8.0 * value) for value in flow_times],
            "plot_t_over_a2_values": sorted_unique([float(row["t_over_a2"]) for row in plot_rows]),
            "plot_t_lat_values": sorted_unique([float(row_t_lat(row)) for row in plot_rows if row_t_lat(row) is not None]),
            "plot_eight_t_over_a2_values": sorted_unique([float(row["eight_t_over_a2"]) for row in plot_rows]),
            "hide_zero_flow_in_plot": bool(args.hide_zero_flow_in_plot),
            "physical_tau_values": [float(value) for value in args.physical_tau],
            "fit_mode": args.fit_mode,
            "fit_excluded_betas": [float(value) for value in fit_excluded_betas],
            "fit_n_points": len(fit_rows),
            "eps_label": eps_label,
            "fixed_epsilon_selection": eps_info,
            "t0_over_a2_by_beta": t0_info,
            "x_axis": x_axis_definition,
            "tau_definition": "tau = t_lat/(t0/a^2)",
            "error_formula": "total_err^2 = stat_err^2 + sys_err^2; sys_err = 0.5*(max_i model_i - min_i model_i)",
            "continuum_error_formula": "When bootstrap refits are available, continuum_err is std_b(c0^(b)) and plotted bands are 16-84 percentiles of refitted curves; otherwise linearized fit covariance is used.",
            "bootstrap_propagation": {
                "enabled": not args.no_continuum_bootstrap,
                "samples_requested": int(args.continuum_bootstrap_samples),
                "seed": int(args.continuum_bootstrap_seed),
                "samples_used": None if bootstrap_matrix is None else int(bootstrap_matrix.shape[1]),
                "row_metadata": bootstrap_metadata,
            },
            "fit_ansatz": {
                "mode": args.fit_mode,
                "formula": fit.get("formula"),
                "fixed_smearing_formula": fit.get("fixed_smearing_formula"),
                "term_names": fit.get("term_names"),
            },
            "continuum_estimates": continuum_estimates,
            "fit": fit,
            "n_points": len(rows),
            "n_points_used_in_fit": len(fit_rows),
            "points": rows,
            "fit_points": fit_rows,
            "skipped": skipped,
            "paths": {
                "table": str(table_path),
                "continuum_limits_table": str(comparison_path),
                "summary": str(summary_path),
                "plot_png": str(continuum_plot_paths.get("png", png_path)),
                "plot_pdf": str(continuum_plot_paths["pdf"]) if "pdf" in continuum_plot_paths else None,
                "plot_html": str(html_path),
                "plot_formats": {fmt: str(fmt_path) for fmt, fmt_path in continuum_plot_paths.items()},
                "bootstrap_npz": str(bootstrap_path) if bootstrap_matrix is not None else None,
            },
            "saved_at": datetime.now(timezone.utc).isoformat(),
        },
    )

    print(f"Saved continuum table: {table_path}")
    print(f"Saved continuum limit comparison: {comparison_path}")
    print(f"Saved continuum summary: {summary_path}")
    print(f"Saved continuum plot: {html_path}")
    if fit_excluded_betas:
        print(
            "Excluded beta values from fit: "
            + ", ".join(f"{value:g}" for value in fit_excluded_betas)
            + f" ({len(rows) - len(fit_rows)} of {len(rows)} point(s))"
        )
    if "pdf" in continuum_plot_paths:
        print(f"Saved continuum plot PDF: {continuum_plot_paths['pdf']}")
    if fit.get("mode") == "combined-risch":
        print(
            "Shared continuum value c00: "
            f"{fit['continuum_value']:.8g} +/- {fit['continuum_err']:.3g} "
            f"(chi2/dof={fit['chi2_dof']})"
        )
        bootstrap_section = fit.get("bootstrap") or {}
        if bootstrap_section:
            print(
                "Bootstrap-propagated c00 statistical error: "
                f"{bootstrap_section.get('continuum_stat_err'):.3g} "
                f"from {bootstrap_section.get('samples_used')} samples"
            )
    else:
        print("Fixed-smearing continuum limits:")
        for estimate in continuum_estimates:
            chi2_text = "nan" if estimate.get("chi2_dof") is None else f"{float(estimate['chi2_dof']):.3g}"
            print(
                f"  {t_lat_label(float(estimate['t_lat']))}: "
                f"c0={float(estimate['continuum_value']):.8g} +/- {float(estimate['continuum_err']):.3g} "
                f"(chi2/dof={chi2_text})"
            )
    if skipped:
        print(f"Skipped {len(skipped)} input/flow item(s); see summary JSON for reasons.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
