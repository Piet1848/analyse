#!/usr/bin/env python3
"""
Search lattice runs by YAML parameters AND calculated results.
"""
from __future__ import annotations
import argparse
import json
import os
import sys
from dataclasses import dataclass
from typing import Any, Dict, Tuple, List, Optional, get_type_hints, get_origin
import concurrent.futures

import numpy as np

from load_input_yaml import load_params, MetropolisParams, GaugeObservableParams
import re

# --- Added tau_int and updated types ---
CALCULATED_FIELDS = {
    "W_R_T": float,
    "W_R_T_err": float,
    "V_R": dict,
    "V_R_err": dict,
    "r0": float,
    "r0_err": float,
    "a": float,
    "a_err": float,
    "epsilon_bar": float,
    "epsilon_bar_err": float,
    "eps_bar": float,
    "eps_bar_err": float,
    "length": float,
    "length_err": float,
    "volume_r0": float,
    "volume_r0_err": float,
    "tau_int": float,
    "block_size": int,
    "analysis_settings": dict,
    "sommer": str,
    "r0_chi": float,
    "r0_chi_err": float,
    "chi": dict,
    "F_chi": dict,
    "creutz_P": dict,
    "creutz_P_err": dict,
    "a_creutz": dict,
    "a_creutz_err": dict,
    "creutz_status": str,
    "gradient_flow": dict,
    "Ehat_clover": dict,
    "Ehat_clover_err": dict,
    "t2E_clover": dict,
    "t2E_clover_err": dict,
    "t0": float,
    "t0_err": float,
}

CALCULATED_FIELD_ALIASES = {
    "eps_bar": "epsilon_bar",
    "eps_bar_err": "epsilon_bar_err",
}

FIELD_ALIASES = {
    "eps1": "epsilon1",
    "eps2": "epsilon2",
}


def build_field_map() -> Dict[str, Tuple[str, Any]]:
    field_map: Dict[str, Tuple[str, Any]] = {}
    metro_types = get_type_hints(MetropolisParams)
    for name, typ in metro_types.items():
        field_map[name] = ("metro", typ)
    gauge_types = get_type_hints(GaugeObservableParams)
    for name, typ in gauge_types.items():
        field_map[name] = ("gauge", typ)
    for name, typ in CALCULATED_FIELDS.items():
        field_map[name] = ("calc", typ)
    return field_map


FIELD_MAP = build_field_map()


def _get_run_evaluation():
    import run_evaluation
    return run_evaluation


def parse_bool(s: str) -> bool:
    sl = s.lower()
    if sl in ("true", "1", "yes", "y", "t"):
        return True
    if sl in ("false", "0", "no", "n", "f"):
        return False
    raise ValueError(f"Cannot parse boolean from '{s}'")


def parse_dynamic_token(tok: str):
    """Parses tokens like V_R5, W_R5_T3, W_R5_T, or W_R_T4. Handles _err suffix."""
    is_err = tok.endswith("_err")
    check_tok = tok[:-4] if is_err else tok

    m = re.fullmatch(r"V_R(\d+)", check_tok)
    if m:
        return ("V_R_err" if is_err else "V_R"), {"R": int(m.group(1))}

    m = re.fullmatch(r"W_R(\d+)_T(\d+)", check_tok)
    if m:
        return ("W_R_T_err" if is_err else "W_R_T"), {"R": int(m.group(1)), "T": int(m.group(2))}

    m = re.fullmatch(r"W_R(\d+)_T", check_tok)
    if m:
        return ("W_R_T_err" if is_err else "W_R_T"), {"R": int(m.group(1))}

    m = re.fullmatch(r"W_R_T(\d+)", check_tok)
    if m:
        return ("W_R_T_err" if is_err else "W_R_T"), {"T": int(m.group(1))}

    m = re.fullmatch(r"creutz_P(\d+)", check_tok)
    if m:
        return ("creutz_P_err" if is_err else "creutz_P"), {"R": int(m.group(1))}

    m = re.fullmatch(r"a_creutz(\d+)", check_tok)
    if m:
        return ("a_creutz_err" if is_err else "a_creutz"), {"R": int(m.group(1))}

    m = re.fullmatch(r"Ehat_clover(.+)", check_tok)
    if m:
        return ("Ehat_clover_err" if is_err else "Ehat_clover"), {"flow_time": m.group(1).replace("p", ".")}

    m = re.fullmatch(r"t2E_clover(.+)", check_tok)
    if m:
        return ("t2E_clover_err" if is_err else "t2E_clover"), {"flow_time": m.group(1).replace("p", ".")}

    return None, None


def get_field_info(name: str):
    """Unified lookup for static and dynamic fields."""
    name = FIELD_ALIASES.get(name, name)
    if name in FIELD_MAP:
        return FIELD_MAP[name]
    base, _ = parse_dynamic_token(name)
    if base:
        return ("calc", Any)
    return None


def convert_value(value_str: str, typ: Any) -> Any:
    origin = get_origin(typ)
    if origin is not None:
        if origin is list or origin is dict or typ is list or typ is dict:
            return value_str
        try:
            return value_str
        except Exception:
            pass
        return value_str

    if typ is int:
        return int(value_str)
    if typ is float:
        return float(value_str)
    if typ is bool:
        return parse_bool(value_str)
    if typ is dict:
        return value_str
    return value_str


def parse_tokens(tokens: list[str]) -> tuple[Dict[str, Any], list[str]]:
    criteria: Dict[str, Any] = {}
    outputs: list[str] = []
    for tok in tokens:
        if "=" in tok:
            name, value_str = tok.split("=", 1)
            name = FIELD_ALIASES.get(name.strip(), name.strip())
            value_str = value_str.strip()
            if name not in FIELD_MAP:
                raise KeyError(f"Unknown parameter '{name}'")
            _, typ = FIELD_MAP[name]
            if FIELD_MAP[name][0] == "calc":
                print(
                    f"Warning: Filtering by calculated field '{name}' triggers calculation.",
                    file=sys.stderr,
                )
            value = convert_value(value_str, typ)
            criteria[name] = value
        else:
            name = tok.strip()
            if not name:
                continue
            info = get_field_info(name)
            if not info:
                raise KeyError(f"Unknown parameter '{name}'")
            outputs.append(name)
    return criteria, outputs


def matches_criteria(
    metro: MetropolisParams,
    gauge: GaugeObservableParams,
    criteria: Dict[str, Any],
) -> bool:
    for name, expected in criteria.items():
        block_name, _ = FIELD_MAP[name]
        if block_name == "calc":
            continue
        obj = metro if block_name == "metro" else gauge
        actual = getattr(obj, name)

        if isinstance(actual, (list, dict)) and isinstance(expected, str):
            if expected not in str(actual):
                return False
            continue

        if actual != expected:
            return False
    return True


def find_matching_runs(root: str, criteria: Dict[str, Any]) -> list[str]:
    matches: list[str] = []
    for dirpath, dirnames, filenames in os.walk(root):
        if "input.yaml" not in filenames:
            continue
        yaml_path = os.path.join(dirpath, "input.yaml")
        try:
            metro, gauge = load_params(yaml_path)
        except Exception as e:
            print(f"Warning: Failed to load parameters from {yaml_path}: {e}", file=sys.stderr)
            continue
        if matches_criteria(metro, gauge, criteria):
            matches.append(dirpath)
    return matches


@dataclass(frozen=True)
class RowSpec:
    row_label: str
    params_path: str
    calc_key: Optional[str]
    is_combined: bool
    n_paths: int = 1


def _get_combined_label(calc_data: Optional[Dict[str, Any]], n_paths: int) -> str:
    if calc_data and "aggregation" in calc_data:
        n_combined = calc_data["aggregation"].get("n_runs_in_group", n_paths)
        return f"combined {n_combined}"
    return f"combined {n_paths}"


def _extract_calculated_value(name: str, calc_data: Optional[Dict[str, Any]]) -> Any:
    if not calc_data or "error" in calc_data:
        return None
    canonical_name = CALCULATED_FIELD_ALIASES.get(name, name)
    if canonical_name in calc_data:
        return calc_data[canonical_name]
    if canonical_name in {"Ehat_clover", "Ehat_clover_err", "t2E_clover", "t2E_clover_err", "t0", "t0_err"}:
        return calc_data.get("gradient_flow", {}).get(canonical_name)

    base, params = parse_dynamic_token(name)
    if not base:
        return None

    if base in ["V_R", "V_R_err", "creutz_P", "creutz_P_err", "a_creutz", "a_creutz_err"]:
        return calc_data.get(base, {}).get(str(params["R"]))

    if base in ["Ehat_clover", "Ehat_clover_err", "t2E_clover", "t2E_clover_err"]:
        flow_dict = calc_data.get("gradient_flow", {}).get(base, {})
        key = str(params["flow_time"])
        try:
            key = f"{float(key):.12g}"
        except (TypeError, ValueError):
            pass
        return flow_dict.get(key)

    if base.startswith("W_R_T"):
        w_dict = calc_data.get(base, {})
        r, t = params.get("R"), params.get("T")
        if r is not None and t is not None:
            return w_dict.get(f"{r},{t}")
        if r is not None:
            sub = {k.split(',')[1]: v for k, v in w_dict.items() if k.startswith(f"{r},")}
            return sub if sub else None
        if t is not None:
            sub = {k.split(',')[0]: v for k, v in w_dict.items() if k.endswith(f",{t}")}
            return sub if sub else None

    return None


def _build_row(
    row_spec: RowSpec,
    outputs: List[str],
    calc_data: Optional[Dict[str, Any]],
) -> Tuple[str, List[Any]]:
    try:
        metro, gauge = load_params(os.path.join(row_spec.params_path, "input.yaml"))
        row_label = row_spec.row_label
        if row_spec.is_combined:
            row_label = _get_combined_label(calc_data, row_spec.n_paths)

        vals = []
        for name in outputs:
            block, _ = get_field_info(name)
            if block == "calc":
                vals.append(_extract_calculated_value(name, calc_data))
            else:
                obj = metro if block == "metro" else gauge
                canonical_name = FIELD_ALIASES.get(name, name)
                vals.append(getattr(obj, canonical_name))
        return row_label, vals
    except Exception:
        row_label = row_spec.row_label
        if row_spec.is_combined:
            row_label = _get_combined_label(calc_data, row_spec.n_paths)
        return row_label, [None] * len(outputs)


def _calculate_task(
    calc_key: str,
    path: str,
    combine_equivalent_runs: bool,
    calc_workers: Optional[int],
    load_workers: int,
    force_calculate: bool = False,
    analysis_options: Optional[Dict[str, Any]] = None,
) -> Tuple[str, Dict[str, Any]]:
    run_evaluation = _get_run_evaluation()
    return calc_key, run_evaluation.get_or_calculate(
        path,
        force_recalc=force_calculate,
        calc_workers=calc_workers,
        load_workers=load_workers,
        combine_equivalent_runs=combine_equivalent_runs,
        analysis_options=analysis_options,
    )


def _print_analysis_settings(
    row_specs: List[RowSpec],
    calc_results: Dict[str, Dict[str, Any]],
) -> None:
    for i, row_spec in enumerate(row_specs):
        calc_data = calc_results.get(row_spec.calc_key) if row_spec.calc_key else None
        row_label = row_spec.row_label
        if row_spec.is_combined:
            row_label = _get_combined_label(calc_data, row_spec.n_paths)

        print(row_label)
        settings = _extract_calculated_value("analysis_settings", calc_data)
        if settings is None:
            print(json.dumps({"error": "analysis settings unavailable"}, indent=2, sort_keys=True))
        else:
            print(json.dumps(settings, indent=2, sort_keys=True))
        if i != len(row_specs) - 1:
            print()


def _extract_row_field(
    row_spec: RowSpec,
    field_name: str,
    calc_data: Optional[Dict[str, Any]],
) -> Any:
    metro, gauge = load_params(os.path.join(row_spec.params_path, "input.yaml"))
    block, _ = get_field_info(field_name)
    if block == "calc":
        return _extract_calculated_value(field_name, calc_data)
    obj = metro if block == "metro" else gauge
    return getattr(obj, FIELD_ALIASES.get(field_name, field_name))


def _coerce_float(value: Any) -> Optional[float]:
    if value is None:
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(out):
        return None
    return out


def _maybe_error_field_name(field_name: str) -> Optional[str]:
    candidates: List[str] = []
    if not field_name.endswith("_err"):
        candidates.append(f"{field_name}_err")
    canonical = CALCULATED_FIELD_ALIASES.get(field_name, field_name)
    if canonical != field_name and not canonical.endswith("_err"):
        candidates.append(f"{canonical}_err")

    for candidate in candidates:
        if get_field_info(candidate):
            return candidate
    return None


def _fit_linear_model(
    x: np.ndarray,
    y: np.ndarray,
    y_err: Optional[np.ndarray] = None,
) -> Tuple[float, float]:
    valid = np.isfinite(x) & np.isfinite(y)
    if y_err is not None:
        valid &= np.isfinite(y_err)

    x_fit = x[valid]
    y_fit = y[valid]
    if x_fit.size < 2 or np.unique(x_fit).size < 2:
        raise ValueError("Need at least two distinct finite x-values for a linear fit.")

    if y_err is not None:
        err_fit = y_err[valid]
        valid_w = err_fit > 0
        if np.count_nonzero(valid_w) >= 2 and np.unique(x_fit[valid_w]).size >= 2:
            p = np.polyfit(x_fit[valid_w], y_fit[valid_w], 1, w=1.0 / err_fit[valid_w])
            return float(p[0]), float(p[1])

    p = np.polyfit(x_fit, y_fit, 1)
    return float(p[0]), float(p[1])


def _bootstrap_linear_model(
    x: np.ndarray,
    y: np.ndarray,
    x_err: Optional[np.ndarray],
    y_err: Optional[np.ndarray],
    n_bootstrap: int,
    seed: int,
) -> Tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    n_points = len(x)
    slopes = np.full(n_bootstrap, np.nan, dtype=float)
    intercepts = np.full(n_bootstrap, np.nan, dtype=float)

    if n_points < 2:
        return slopes, intercepts

    for i in range(n_bootstrap):
        sample_idx = rng.integers(0, n_points, size=n_points)
        x_sample = np.asarray(x[sample_idx], dtype=float)
        y_sample = np.asarray(y[sample_idx], dtype=float)
        x_sigma = None if x_err is None else np.asarray(x_err[sample_idx], dtype=float)
        y_sigma = None if y_err is None else np.asarray(y_err[sample_idx], dtype=float)

        if x_sigma is not None:
            valid_x_sigma = np.isfinite(x_sigma) & (x_sigma > 0)
            if np.any(valid_x_sigma):
                x_sample[valid_x_sigma] = rng.normal(x_sample[valid_x_sigma], x_sigma[valid_x_sigma])
        if y_sigma is not None:
            valid_y_sigma = np.isfinite(y_sigma) & (y_sigma > 0)
            if np.any(valid_y_sigma):
                y_sample[valid_y_sigma] = rng.normal(y_sample[valid_y_sigma], y_sigma[valid_y_sigma])

        try:
            slope, intercept = _fit_linear_model(x_sample, y_sample, y_sigma)
        except (ValueError, np.linalg.LinAlgError):
            continue

        slopes[i] = slope
        intercepts[i] = intercept

    return slopes, intercepts


def _summarize_bootstrap(values: np.ndarray) -> Tuple[Optional[float], Optional[float], Optional[float], Optional[float]]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return None, None, None, None
    return (
        float(np.mean(finite)),
        float(np.std(finite)),
        float(np.percentile(finite, 16)),
        float(np.percentile(finite, 84)),
    )


def _print_linear_fit_report(
    row_specs: List[RowSpec],
    calc_results: Dict[str, Dict[str, Any]],
    x_name: str,
    y_name: str,
    exclude_origin: bool,
    n_bootstrap: int,
    seed: int,
    predict_x: List[float],
    predict_y: List[float],
) -> None:
    x_err_name = _maybe_error_field_name(x_name)
    y_err_name = _maybe_error_field_name(y_name)

    fit_rows: List[Tuple[str, float, float, Optional[float], Optional[float]]] = []
    for row_spec in row_specs:
        calc_data = calc_results.get(row_spec.calc_key) if row_spec.calc_key else None
        row_label = row_spec.row_label
        if row_spec.is_combined:
            row_label = _get_combined_label(calc_data, row_spec.n_paths)

        try:
            x_val = _coerce_float(_extract_row_field(row_spec, x_name, calc_data))
            y_val = _coerce_float(_extract_row_field(row_spec, y_name, calc_data))
            x_err = _coerce_float(_extract_row_field(row_spec, x_err_name, calc_data)) if x_err_name else None
            y_err = _coerce_float(_extract_row_field(row_spec, y_err_name, calc_data)) if y_err_name else None
        except Exception:
            continue

        if x_val is None or y_val is None:
            continue
        if exclude_origin and np.isclose(x_val, 0.0) and np.isclose(y_val, 0.0):
            continue
        fit_rows.append((row_label, x_val, y_val, x_err, y_err))

    if len(fit_rows) < 2:
        raise ValueError("Not enough finite points remain for the linear fit.")

    x = np.asarray([row[1] for row in fit_rows], dtype=float)
    y = np.asarray([row[2] for row in fit_rows], dtype=float)
    x_err = np.asarray(
        [np.nan if row[3] is None else row[3] for row in fit_rows],
        dtype=float,
    ) if x_err_name else None
    y_err = np.asarray(
        [np.nan if row[4] is None else row[4] for row in fit_rows],
        dtype=float,
    ) if y_err_name else None

    slope, intercept = _fit_linear_model(x, y, y_err)
    slopes_boot, intercepts_boot = _bootstrap_linear_model(
        x,
        y,
        x_err,
        y_err,
        n_bootstrap=n_bootstrap,
        seed=seed,
    )
    slope_mean, slope_std, slope_p16, slope_p84 = _summarize_bootstrap(slopes_boot)
    int_mean, int_std, int_p16, int_p84 = _summarize_bootstrap(intercepts_boot)

    print()
    print(f"Linear fit ({y_name} = m * {x_name} + b)")
    print(f"  points used: {len(fit_rows)}")
    print(f"  excluded (0,0): {'yes' if exclude_origin else 'no'}")
    print(f"  x error field: {x_err_name if x_err_name else 'none'}")
    print(f"  y error field: {y_err_name if y_err_name else 'none'}")
    print(f"  bootstrap replicas: {n_bootstrap}")
    if slope_std is not None:
        print(
            f"  slope m: {slope:.12g}  (bootstrap mean {slope_mean:.12g}, std {slope_std:.6g}, 16-84% [{slope_p16:.12g}, {slope_p84:.12g}])"
        )
    else:
        print(f"  slope m: {slope:.12g}")
    if int_std is not None:
        print(
            f"  intercept b: {intercept:.12g}  (bootstrap mean {int_mean:.12g}, std {int_std:.6g}, 16-84% [{int_p16:.12g}, {int_p84:.12g}])"
        )
    else:
        print(f"  intercept b: {intercept:.12g}")

    if predict_y:
        print()
        print(f"Predicted {y_name} from {x_name}")
        print(f"{x_name:<16}  {y_name:<18}  {y_name}_boot_err")
        for x_target in predict_y:
            y_central = slope * x_target + intercept
            y_boot = slopes_boot * x_target + intercepts_boot
            _, y_std, _, _ = _summarize_bootstrap(y_boot)
            y_std_str = f"{y_std:.6g}" if y_std is not None else "N/A"
            print(f"{x_target:<16.12g}  {y_central:<18.12g}  {y_std_str}")

    if predict_x:
        print()
        print(f"Predicted {x_name} from {y_name}")
        print(f"{y_name:<16}  {x_name:<18}  {x_name}_boot_err")
        slope_mask = np.isfinite(slopes_boot) & (~np.isclose(slopes_boot, 0.0))
        slopes_inv = slopes_boot[slope_mask]
        intercepts_inv = intercepts_boot[slope_mask]
        for y_target in predict_x:
            x_central = np.nan if np.isclose(slope, 0.0) else (y_target - intercept) / slope
            if slopes_inv.size > 0:
                x_boot = (y_target - intercepts_inv) / slopes_inv
            else:
                x_boot = np.asarray([], dtype=float)
            _, x_std, _, _ = _summarize_bootstrap(x_boot)
            x_val_str = f"{x_central:.12g}" if np.isfinite(x_central) else "N/A"
            x_std_str = f"{x_std:.6g}" if x_std is not None else "N/A"
            print(f"{y_target:<16.12g}  {x_val_str:<18}  {x_std_str}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Search data runs and calculate missing values.")
    parser.add_argument("root", help="Root data directory")
    parser.add_argument("tokens", nargs="*", help="Filters (NAME=VAL) and outputs (NAME)")
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of runs to calculate in parallel (higher values use more RAM)",
    )
    parser.add_argument(
        "--calc-workers",
        type=int,
        default=None,
        help="Threads per run for W(R,T) and V(R) calculations",
    )
    parser.add_argument(
        "--load-workers",
        type=int,
        default=1,
        help="Files to load in parallel while combining equivalent runs",
    )
    parser.add_argument(
        "--force-calculate",
        "--force_calculate",
        dest="force_calculate",
        action="store_true",
        help="Recalculate results instead of loading cached analysis data",
    )
    parser.add_argument(
        "--scope",
        choices=["both", "single", "combined"],
        default="both",
        help="Which rows to emit: individual runs, combined groups, or both",
    )
    parser.add_argument(
        "--show-analysis-settings",
        action="store_true",
        help="Print the key analysis settings used for each matching row",
    )
    parser.add_argument(
        "--r0-t-min",
        type=int,
        default=None,
        help="Override the t_min used for the r0 fit",
    )
    parser.add_argument(
        "--fit-linear",
        nargs=2,
        metavar=("X", "Y"),
        help="Fit a linear model Y = m * X + b using the emitted rows.",
    )
    parser.add_argument(
        "--fit-bootstrap",
        type=int,
        default=2000,
        help="Bootstrap replicas for --fit-linear (default: 2000).",
    )
    parser.add_argument(
        "--fit-seed",
        type=int,
        default=42,
        help="Random seed for --fit-linear bootstrap (default: 42).",
    )
    parser.add_argument(
        "--fit-exclude-origin",
        action="store_true",
        help="Exclude rows where both fit coordinates are exactly zero.",
    )
    parser.add_argument(
        "--predict-y",
        nargs="*",
        type=float,
        default=[],
        metavar="X",
        help="After --fit-linear, evaluate the fitted Y at one or more X values.",
    )
    parser.add_argument(
        "--predict-x",
        nargs="*",
        type=float,
        default=[],
        metavar="Y",
        help="After --fit-linear, invert the fitted relation and estimate X for one or more Y values.",
    )
    args = parser.parse_args()

    if not args.tokens and not args.show_analysis_settings:
        print("Usage:")
        print("  python search_data.py ROOT [TOKENS ...] [--workers N] [--calc-workers N] [--load-workers N]")
        print("  python search_data.py ROOT --show-analysis-settings")
        print()
        print("Examples:")
        print("  python search_data.py ../data beta=2.4 lattice_size")
        print("  python search_data.py ../data W_R_T --workers 4 --calc-workers 2 --load-workers 4")
        print("  python search_data.py ../data beta=2.4 --show-analysis-settings")
        print()
        print("Options:")
        print("  --workers N         Number of runs to calculate in parallel (default: 1)")
        print("  --calc-workers N    Threads per run for W(R,T) and V(R) calculations (default: auto)")
        print("  --load-workers N    Files to load in parallel while combining equivalent runs (default: 1)")
        print("  --force-calculate   Recalculate results instead of loading cached analysis data")
        print("  --scope {both,single,combined}  Which rows to emit (default: both)")
        print("  --show-analysis-settings  Print the key analysis settings JSON per row")
        print("  --r0-t-min N        Override the t_min used for the r0 fit")
        print("  --fit-linear X Y    Fit Y = m * X + b over the emitted rows")
        print("  --fit-bootstrap N   Bootstrap replicas for --fit-linear (default: 2000)")
        print("  --fit-seed N        Random seed for the linear-fit bootstrap")
        print("  --fit-exclude-origin  Exclude rows where X=0 and Y=0 from the fit")
        print("  --predict-y X ...   Evaluate fitted Y values for one or more X inputs")
        print("  --predict-x Y ...   Invert fitted relation for one or more Y inputs")
        print()
        print("Available Parameters:")
        for k, v in FIELD_MAP.items():
            if v[0] != "calc":
                print(f"  {k}")
        print("\nCalculated Fields:")
        for k, v in FIELD_MAP.items():
            if v[0] == "calc":
                print(f"  {k}")
        return

    try:
        criteria, outputs = parse_tokens(args.tokens)
    except Exception as e:
        sys.exit(f"Error: {e}")

    if args.fit_linear:
        fit_x, fit_y = args.fit_linear
        if not get_field_info(fit_x):
            sys.exit(f"Error: Unknown fit field '{fit_x}'")
        if not get_field_info(fit_y):
            sys.exit(f"Error: Unknown fit field '{fit_y}'")
        if args.fit_bootstrap < 1:
            sys.exit("Error: --fit-bootstrap must be at least 1")
    elif args.predict_x or args.predict_y:
        sys.exit("Error: --predict-x and --predict-y require --fit-linear")

    if args.show_analysis_settings and "analysis_settings" not in outputs:
        outputs.append("analysis_settings")

    matches = find_matching_runs(args.root, criteria)
    if not matches:
        print("No matching runs found.")
        return

    from collections import defaultdict

    groups = defaultdict(list)
    run_evaluation = _get_run_evaluation()
    for path in matches:
        group_key = run_evaluation._group_key_for_run(path)
        if group_key is None:
            groups[path].append(path)
        else:
            groups[group_key].append(path)

    row_specs: List[RowSpec] = []
    calc_targets: Dict[str, Tuple[str, bool]] = {}
    for group_key, paths in groups.items():
        representative_path = paths[0]
        include_single = args.scope in ("both", "single")
        include_combined = args.scope in ("both", "combined") and len(paths) > 1

        if include_single:
            for path in paths:
                calc_key = f"single:{path}"
                row_specs.append(
                    RowSpec(
                        row_label=path,
                        params_path=path,
                        calc_key=calc_key,
                        is_combined=False,
                    )
                )
                calc_targets[calc_key] = (path, False)

        if include_combined:
            calc_key = f"combined:{representative_path}"
            row_specs.append(
                RowSpec(
                    row_label=f"combined {len(paths)}",
                    params_path=representative_path,
                    calc_key=calc_key,
                    is_combined=True,
                    n_paths=len(paths),
                )
            )
            calc_targets[calc_key] = (representative_path, True)

    if not row_specs:
        print(f"No rows found for scope '{args.scope}'.")
        return

    n_combined_groups = sum(1 for p in groups.values() if len(p) > 1)
    needs_calc = args.force_calculate or any(get_field_info(out)[0] == "calc" for out in outputs)
    calc_results: Dict[str, Dict[str, Any]] = {}
    analysis_options: Dict[str, Any] = {}
    if args.r0_t_min is not None:
        analysis_options["r0_t_min"] = args.r0_t_min

    if needs_calc:
        print(
            f"Processing {len(row_specs)} rows ({len(matches)} individual runs, {n_combined_groups} combined groups) "
            f"with {len(calc_targets)} unique calculation task(s), "
            f"{args.workers} run worker(s), "
            f"{args.calc_workers if args.calc_workers is not None else 'auto'} calc worker(s) per run, "
            f"and {args.load_workers} load worker(s)"
            f"{' with cache reuse disabled' if args.force_calculate else ''}...",
            file=sys.stderr,
        )
        if args.workers <= 1:
            for i, (calc_key, (path, combine_equivalent_runs)) in enumerate(calc_targets.items(), start=1):
                result_key, calc_data = _calculate_task(
                    calc_key,
                    path,
                    combine_equivalent_runs,
                    args.calc_workers,
                    args.load_workers,
                    args.force_calculate,
                    analysis_options=analysis_options,
                )
                calc_results[result_key] = calc_data
                print(f"Calculated {i}/{len(calc_targets)}", file=sys.stderr, end="\r")
        else:
            with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as executor:
                futures = {
                    executor.submit(
                        _calculate_task,
                        calc_key,
                        path,
                        combine_equivalent_runs,
                        args.calc_workers,
                        args.load_workers,
                        args.force_calculate,
                        analysis_options,
                    ): calc_key
                    for calc_key, (path, combine_equivalent_runs) in calc_targets.items()
                }
                for i, future in enumerate(concurrent.futures.as_completed(futures), start=1):
                    result_key, calc_data = future.result()
                    calc_results[result_key] = calc_data
                    print(f"Calculated {i}/{len(calc_targets)}", file=sys.stderr, end="\r")
        print(" " * 60, file=sys.stderr, end="\r")
    else:
        print(
            f"Processing {len(row_specs)} rows ({len(matches)} individual runs, {n_combined_groups} combined groups) "
            f"without calculated fields...",
            file=sys.stderr,
        )

    if outputs == ["analysis_settings"]:
        _print_analysis_settings(row_specs, calc_results)
        return

    rows = []
    for row_spec in row_specs:
        calc_data = calc_results.get(row_spec.calc_key) if row_spec.calc_key else None
        rows.append(_build_row(row_spec, outputs, calc_data))

    def sort_key(row):
        key_parts = []
        for val in row[1]:
            if val is None:
                key_parts.append((0, ""))
            elif isinstance(val, (dict, list)):
                key_parts.append((1, str(val)))
            else:
                key_parts.append((1, val))
        return tuple(key_parts)

    rows.sort(key=sort_key)

    headers = ["path"] + outputs
    all_rows = []
    for row_label, vals in rows:
        out_strs = [str(v) if v is not None else "N/A" for v in vals]
        all_rows.append([row_label] + out_strs)

    col_widths = [len(h) for h in headers]
    for row in all_rows:
        for i, val in enumerate(row):
            col_widths[i] = max(col_widths[i], len(val))

    fmt = "  ".join(f"{{:<{w}}}" for w in col_widths)

    print(fmt.format(*headers))
    for row in all_rows:
        print(fmt.format(*row))

    if args.fit_linear:
        fit_x, fit_y = args.fit_linear
        _print_linear_fit_report(
            row_specs=row_specs,
            calc_results=calc_results,
            x_name=fit_x,
            y_name=fit_y,
            exclude_origin=args.fit_exclude_origin,
            n_bootstrap=args.fit_bootstrap,
            seed=args.fit_seed,
            predict_x=args.predict_x,
            predict_y=args.predict_y,
        )


if __name__ == "__main__":
    main()
