#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import html
import os
import shlex
import shutil
import subprocess
import sys
from dataclasses import asdict
from datetime import datetime, timezone
from math import ceil
from pathlib import Path
from typing import Any, Callable, get_origin, get_type_hints

import numpy as np

import data_organizer
import run_evaluation
from calculator import Calculator, fit_r0_from_potential_data
from finalized_analysis_helpers import (
    build_bootstrap_block_size_scan,
    build_effective_mass_scan,
    build_thermalization_preview,
    build_locked_r0_scan,
    build_v_r_scan,
    build_wrt_scan,
    enumerate_t_windows,
    load_json,
    recommend_bootstrap_block_size,
    save_bootstrap_block_size_plot,
    save_cornell_plot,
    save_creutz_diagonal_plot,
    save_creutz_plot,
    save_effective_mass_plot,
    save_gradient_flow_plot,
    save_json,
    save_thermalization_plot,
    save_r0_stability_plot,
    summarize_bootstrap,
    window_label,
)
from load_input_yaml import (
    GaugeObservableParams,
    GradientFlowParams,
    MetropolisParams,
    load_gradient_flow_params,
    load_params,
)


SCHEMA_VERSION = 1
IGNORED_METRO_FIELDS = {"seed", "nSweep"}
REQUIRED_METRO_FIELDS = {"beta", "L0", "epsilon1"}
ALLOWED_GRADIENT_FLOW_FIELDS = {"dt", "t_values"}
GROUP_FIELD_ALIASES = {
    "eps1": "epsilon1",
    "dt": "gradient_flow.dt",
    "t_values": "gradient_flow.t_values",
    "gf.dt": "gradient_flow.dt",
    "gradient.dt": "gradient_flow.dt",
    "gradient_flow_dt": "gradient_flow.dt",
    "gradient_flow_t_values": "gradient_flow.t_values",
}
FILTER_FIELD_ALIASES = GROUP_FIELD_ALIASES
DEFAULT_BLOCK_SIZE_SCAN_MIN = 1
DEFAULT_BLOCK_SIZE_SCAN_MAX = None
DEFAULT_BLOCK_SIZE_SCAN_MAX_FRACTION = 16
DEFAULT_BLOCK_SIZE_SCAN_STEP = 1


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def format_token(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:g}"
    return str(value)


def format_result_value(value: Any, error: Any = None) -> str:
    if value is None:
        return "n/a"
    try:
        numeric_value = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not np.isfinite(numeric_value):
        return str(numeric_value)

    value_text = f"{numeric_value:.6g}"
    if error is None:
        return value_text
    try:
        numeric_error = float(error)
    except (TypeError, ValueError):
        return value_text
    if not np.isfinite(numeric_error):
        return value_text
    return f"{value_text} +/- {numeric_error:.3g}"


def uncertainty_decimal_places(error: float) -> int | None:
    if not np.isfinite(error) or error <= 0:
        return None

    exponent = int(np.floor(np.log10(abs(error))))
    leading_digit = int(abs(error) / (10 ** exponent))
    significant_digits = 2 if leading_digit == 1 else 1
    return significant_digits - 1 - exponent


def format_decimal_places(value: float, decimal_places: int) -> str:
    rounded = round(value, decimal_places)
    if rounded == 0:
        rounded = 0.0
    if decimal_places > 0:
        return f"{rounded:.{decimal_places}f}"
    return f"{rounded:.0f}"


def format_result_value_rounded_to_error(value: Any, error: Any) -> str:
    if value is None:
        return "n/a"
    try:
        numeric_value = float(value)
        numeric_error = float(error)
    except (TypeError, ValueError):
        return format_result_value(value, error)
    if not np.isfinite(numeric_value) or not np.isfinite(numeric_error):
        return format_result_value(value, error)

    decimal_places = uncertainty_decimal_places(numeric_error)
    if decimal_places is None:
        return format_result_value(value, error)

    return (
        f"{format_decimal_places(numeric_value, decimal_places)} +/- "
        f"{format_decimal_places(numeric_error, decimal_places)}"
    )


def format_compact_float(value: Any) -> str:
    if value is None:
        return "n/a"
    try:
        numeric_value = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not np.isfinite(numeric_value):
        return str(numeric_value)
    if numeric_value == 0.0:
        return "0.0"
    return f"{numeric_value:.6g}"


def build_filter_field_map() -> dict[str, tuple[str, Any]]:
    field_map: dict[str, tuple[str, Any]] = {}
    for name, typ in get_type_hints(MetropolisParams).items():
        field_map[name] = ("metro", typ)
    for name, typ in get_type_hints(GaugeObservableParams).items():
        field_map[name] = ("gauge", typ)
    for name, typ in get_type_hints(GradientFlowParams).items():
        prefixed_name = f"gradient_flow.{name}"
        field_map[prefixed_name] = ("gradient_flow", typ)
        if name not in field_map:
            field_map[name] = ("gradient_flow", typ)
    return field_map


FILTER_FIELD_MAP = build_filter_field_map()


def parse_bool(value: str) -> bool:
    lowered = value.lower()
    if lowered in {"true", "1", "yes", "y", "t"}:
        return True
    if lowered in {"false", "0", "no", "n", "f"}:
        return False
    raise ValueError(f"Cannot parse boolean from {value!r}")


def convert_filter_value(value: str, typ: Any) -> Any:
    origin = get_origin(typ)
    if origin in {list, dict} or typ in {list, dict}:
        return value
    if typ is bool:
        return parse_bool(value)
    if typ is int:
        return int(value)
    if typ is float:
        return float(value)
    return value


def parse_filter_tokens(tokens: list[str] | None) -> dict[str, Any]:
    criteria: dict[str, Any] = {}
    for token in tokens or []:
        if "=" not in token:
            raise ValueError(f"Filter token must use NAME=VALUE syntax: {token!r}")
        raw_name, raw_value = token.split("=", 1)
        name = FILTER_FIELD_ALIASES.get(raw_name.strip(), raw_name.strip())
        if name not in FILTER_FIELD_MAP:
            raise ValueError(f"Unknown filter field: {raw_name.strip()!r}")
        _, typ = FILTER_FIELD_MAP[name]
        criteria[name] = convert_filter_value(raw_value.strip(), typ)
    return criteria


def run_matches_filter(
    metro: MetropolisParams,
    gauge: GaugeObservableParams,
    gradient_flow: GradientFlowParams | None,
    criteria: dict[str, Any],
) -> bool:
    for name, expected in criteria.items():
        block_name, _ = FILTER_FIELD_MAP[name]
        if block_name == "metro":
            obj = metro
        elif block_name == "gauge":
            obj = gauge
        else:
            obj = gradient_flow
        if obj is None:
            return False
        attr_name = name.split(".", 1)[1] if block_name == "gradient_flow" and "." in name else name
        actual = getattr(obj, attr_name)
        if isinstance(actual, (list, dict)) and isinstance(expected, str):
            if expected not in str(actual):
                return False
            continue
        if actual != expected:
            return False
    return True


def normalize_thermalization_steps_by_run(
    run_dirs: list[str],
    thermalization_steps_by_run: dict[str, int],
) -> dict[str, int]:
    normalized = {os.path.abspath(path): int(value) for path, value in thermalization_steps_by_run.items()}
    ordered: dict[str, int] = {}
    for run_dir in [os.path.abspath(path) for path in run_dirs]:
        if run_dir not in normalized:
            raise ValueError(f"Missing thermalization cut for run: {run_dir}")
        if int(normalized[run_dir]) < 0:
            raise ValueError(f"Thermalization cut must be non-negative for run: {run_dir}")
        ordered[run_dir] = int(normalized[run_dir])
    return ordered


def common_thermalization_step(thermalization_steps_by_run: dict[str, int]) -> int | None:
    unique = sorted({int(value) for value in thermalization_steps_by_run.values()})
    return unique[0] if len(unique) == 1 else None


def wilson_flow_time_pair_filter(
    requested_flow_time: float | None,
) -> Callable[[float | None, int, int], bool]:
    normalized_requested = data_organizer._normalize_flow_time(requested_flow_time)

    def pair_filter(flow_time: float | None, _r_val: int, _t_val: int) -> bool:
        return data_organizer._normalize_flow_time(flow_time) == normalized_requested

    return pair_filter


def wilson_flow_time_filter_label(requested_flow_time: float | None) -> str:
    normalized_requested = data_organizer._normalize_flow_time(requested_flow_time)
    if normalized_requested is None:
        return "flow_time=unflowed"
    return f"flow_time={float(normalized_requested):g}"


def analysis_hash(
    run_dirs: list[str],
    thermalization_steps_by_run: dict[str, int],
    wilson_flow_time: float | None = None,
) -> str:
    normalized_cuts = normalize_thermalization_steps_by_run(run_dirs, thermalization_steps_by_run)
    payload = {
        "run_dirs": sorted(os.path.abspath(path) for path in run_dirs),
        "thermalization_steps_by_run": normalized_cuts,
    }
    if wilson_flow_time is not None:
        payload["wilson_flow_time"] = float(wilson_flow_time)
    digest = hashlib.md5(repr(payload).encode("utf-8")).hexdigest()
    return digest[:12]


def preview_hash(run_dirs: list[str]) -> str:
    payload = {
        "run_dirs": sorted(os.path.abspath(path) for path in run_dirs),
    }
    digest = hashlib.md5(repr(payload).encode("utf-8")).hexdigest()
    return digest[:12]


def thermalization_preview_filename(run_dir: str) -> str:
    run_name = Path(run_dir).name or "run"
    safe_name = "".join(ch if ch.isalnum() or ch in {"-", "_"} else "_" for ch in run_name)
    digest = hashlib.md5(os.path.abspath(run_dir).encode("utf-8")).hexdigest()[:12]
    return f"thermalization_{safe_name}_{digest}.html"


def parse_bootstrap_block_sizes(args: argparse.Namespace) -> list[int] | None:
    if args.block_sizes:
        values: list[int] = []
        for chunk in args.block_sizes.split(","):
            stripped = chunk.strip()
            if not stripped:
                continue
            values.append(int(stripped))
    else:
        min_block = int(args.min_block_size)
        step = int(args.block_step)
        if min_block < 1:
            raise ValueError("--min-block-size must be at least 1")
        if step < 1:
            raise ValueError("--block-step must be at least 1")
        if args.max_block_size is None:
            return None
        max_block = int(args.max_block_size)
        if max_block < min_block:
            raise ValueError("--max-block-size must be greater than or equal to --min-block-size")
        values = list(range(min_block, max_block + 1, step))
        if values and values[-1] != max_block:
            values.append(max_block)

    block_sizes = sorted({int(value) for value in values if int(value) >= 1})
    if not block_sizes:
        raise ValueError("No valid bootstrap block sizes selected.")
    return block_sizes


def open_html_plot(path: Path) -> bool:
    resolved = path.resolve()
    uri = resolved.as_uri()

    if os.name == "nt":
        try:
            os.startfile(str(resolved))  # type: ignore[attr-defined]
            return True
        except Exception:
            return False

    candidates: list[list[str]] = []
    browser_env = os.environ.get("BROWSER", "").strip()
    if browser_env:
        for entry in browser_env.split(os.pathsep):
            tokens = shlex.split(entry)
            if not tokens:
                continue
            if any("%s" in token for token in tokens):
                candidates.append([token.replace("%s", uri) for token in tokens])
            else:
                candidates.append([*tokens, uri])

    if sys.platform == "darwin":
        opener = shutil.which("open")
        if opener:
            candidates.append([opener, str(resolved)])
    else:
        opener = shutil.which("xdg-open")
        if opener:
            candidates.append([opener, str(resolved)])
        gio = shutil.which("gio")
        if gio:
            candidates.append([gio, "open", uri])

    seen: set[tuple[str, ...]] = set()
    for command in candidates:
        key = tuple(command)
        if key in seen:
            continue
        seen.add(key)
        try:
            process = subprocess.Popen(
                command,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                start_new_session=True,
            )
            try:
                returncode = process.wait(timeout=1.0)
            except subprocess.TimeoutExpired:
                return True
            if returncode == 0:
                return True
        except Exception:
            continue

    return False


def compare_dataclass_dicts(reference: dict[str, Any], current: dict[str, Any]) -> dict[str, tuple[Any, Any]]:
    diff: dict[str, tuple[Any, Any]] = {}
    all_keys = sorted(set(reference) | set(current))
    for key in all_keys:
        if reference.get(key) != current.get(key):
            diff[key] = (reference.get(key), current.get(key))
    return diff


def canonical_group_field(name: str) -> str:
    return GROUP_FIELD_ALIASES.get(name, name)


def freeze_group_value(value: Any) -> Any:
    if isinstance(value, dict):
        return tuple((key, freeze_group_value(val)) for key, val in sorted(value.items()))
    if isinstance(value, list):
        return tuple(freeze_group_value(item) for item in value)
    return value


def _jsonable_group_value(value: Any) -> Any:
    if isinstance(value, tuple):
        return [_jsonable_group_value(item) for item in value]
    if isinstance(value, list):
        return [_jsonable_group_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _jsonable_group_value(val) for key, val in value.items()}
    return value


def _unique_values_for_runs(run_values: dict[str, Any]) -> list[Any]:
    frozen_seen: dict[Any, Any] = {}
    for value in run_values.values():
        frozen = freeze_group_value(value)
        frozen_seen.setdefault(frozen, value)
    return [
        _jsonable_group_value(value)
        for _, value in sorted(frozen_seen.items(), key=lambda item: repr(item[0]))
    ]


def _varying_field_summary(run_values: dict[str, dict[str, Any]], field_names: set[str]) -> dict[str, Any]:
    summary: dict[str, Any] = {}
    for field in sorted(field_names):
        values_by_run = {run_dir: values.get(field) for run_dir, values in run_values.items()}
        unique_values = _unique_values_for_runs(values_by_run)
        if len(unique_values) <= 1:
            continue
        summary[field] = {
            "unique_values": unique_values,
            "values_by_run": {
                run_dir: _jsonable_group_value(value)
                for run_dir, value in values_by_run.items()
            },
        }
    return summary


def _common_field_summary(run_values: dict[str, dict[str, Any]], field_names: set[str]) -> dict[str, Any]:
    summary: dict[str, Any] = {}
    for field in sorted(field_names):
        values_by_run = {run_dir: values.get(field) for run_dir, values in run_values.items()}
        unique_values = _unique_values_for_runs(values_by_run)
        if len(unique_values) == 1:
            summary[field] = unique_values[0]
    return summary


def discover_run_directories(root: str) -> list[str]:
    root_path = Path(root).expanduser().resolve()
    if not root_path.is_dir():
        raise ValueError(f"Run root not found: {root_path}")

    run_dirs = sorted(
        {
            str(path.parent.resolve())
            for path in root_path.rglob("input.yaml")
            if path.is_file()
        }
    )
    if not run_dirs:
        raise ValueError(f"No run directories with input.yaml found under {root_path}")
    return run_dirs


def discover_run_directories_from_roots(roots: list[str]) -> list[str]:
    run_dirs: set[str] = set()
    errors: list[str] = []
    for root in roots:
        try:
            run_dirs.update(discover_run_directories(root))
        except ValueError as exc:
            errors.append(str(exc))
    if errors and not run_dirs:
        raise ValueError("; ".join(errors))
    return sorted(run_dirs)


def resolve_selected_run_directories(paths: list[str]) -> tuple[list[str], bool]:
    run_dirs: set[str] = set()
    expanded_roots = False
    missing: list[str] = []
    errors: list[str] = []

    for raw_path in paths:
        path = Path(raw_path).expanduser().resolve()
        if not path.is_dir():
            missing.append(str(path))
            continue
        if (path / "input.yaml").is_file():
            run_dirs.add(str(path))
            continue
        expanded_roots = True
        try:
            run_dirs.update(discover_run_directories(str(path)))
        except ValueError as exc:
            errors.append(str(exc))

    if missing:
        raise ValueError(f"Run directories not found: {missing}")
    if errors and not run_dirs:
        raise ValueError("; ".join(errors))
    return sorted(run_dirs), expanded_roots


def filter_run_directories(run_dirs: list[str], criteria: dict[str, Any]) -> list[str]:
    if not criteria:
        return sorted(os.path.abspath(path) for path in run_dirs)

    matches: list[str] = []
    for run_dir in sorted(os.path.abspath(path) for path in run_dirs):
        yaml_path = Path(run_dir) / "input.yaml"
        if not yaml_path.exists():
            continue
        metro, gauge = load_params(str(yaml_path))
        gradient_flow = load_gradient_flow_params(str(yaml_path))
        if run_matches_filter(metro, gauge, gradient_flow, criteria):
            matches.append(run_dir)
    return matches


def exclude_run_directories(run_dirs: list[str], patterns: list[str] | None) -> list[str]:
    if not patterns:
        return sorted(os.path.abspath(path) for path in run_dirs)

    excluded: list[str] = []
    for run_dir in sorted(os.path.abspath(path) for path in run_dirs):
        path = Path(run_dir)
        if any(pattern in run_dir or path.match(pattern) for pattern in patterns):
            continue
        excluded.append(run_dir)
    return excluded


def build_group_signature(
    metro: MetropolisParams,
    gauge: GaugeObservableParams,
    gradient_flow: GradientFlowParams | None = None,
    group_by: list[str] | None = None,
    require_fixed_dt: bool = False,
) -> tuple[tuple[str, Any], ...]:
    all_metro_fields = {k: v for k, v in asdict(metro).items() if k not in IGNORED_METRO_FIELDS}
    metro_fields = {k: v for k, v in all_metro_fields.items() if k in REQUIRED_METRO_FIELDS}
    gauge_fields = asdict(gauge)
    ignored_gradient_flow_fields = set(ALLOWED_GRADIENT_FLOW_FIELDS)
    if require_fixed_dt:
        ignored_gradient_flow_fields.discard("dt")
    all_gradient_flow_fields = (
        {f"gradient_flow.{key}": value for key, value in asdict(gradient_flow).items()}
        if gradient_flow is not None
        else {}
    )

    if group_by:
        combined = {**all_metro_fields, **gauge_fields, **all_gradient_flow_fields}
        requested = [canonical_group_field(name) for name in group_by]
        missing = [name for name in requested if name not in combined]
        if missing:
            raise ValueError(f"Unknown grouping field(s): {', '.join(missing)}")
        return tuple((name, freeze_group_value(combined[name])) for name in requested)

    gradient_flow_fields = {}
    if require_fixed_dt and "gradient_flow.dt" in all_gradient_flow_fields:
        gradient_flow_fields["gradient_flow.dt"] = all_gradient_flow_fields["gradient_flow.dt"]
    combined = {**metro_fields, **gradient_flow_fields}
    return tuple((name, freeze_group_value(value)) for name, value in sorted(combined.items()))


def describe_group_signature(signature: tuple[tuple[str, Any], ...]) -> str:
    return ", ".join(f"{name}={format_token(value)}" for name, value in signature)


def group_run_directories(
    run_dirs: list[str],
    group_by: list[str] | None = None,
    require_fixed_dt: bool = False,
) -> list[tuple[tuple[tuple[str, Any], ...], list[str]]]:
    grouped: dict[tuple[tuple[str, Any], ...], list[str]] = {}

    for run_dir in sorted(os.path.abspath(path) for path in run_dirs):
        yaml_path = Path(run_dir) / "input.yaml"
        if not yaml_path.exists():
            raise ValueError(f"Missing input.yaml in {run_dir}")
        metro, gauge = load_params(str(yaml_path))
        gradient_flow = load_gradient_flow_params(str(yaml_path))
        signature = build_group_signature(
            metro,
            gauge,
            gradient_flow,
            group_by=group_by,
            require_fixed_dt=require_fixed_dt,
        )
        grouped.setdefault(signature, []).append(run_dir)

    return sorted(grouped.items(), key=lambda item: (describe_group_signature(item[0]), item[1]))


def validate_run_directories(
    run_dirs: list[str],
    *,
    require_fixed_dt: bool = False,
) -> tuple[MetropolisParams, GaugeObservableParams, dict[str, Any]]:
    if not run_dirs:
        raise ValueError("At least one run directory is required.")

    normalized = [os.path.abspath(path) for path in run_dirs]
    missing = [path for path in normalized if not os.path.isdir(path)]
    if missing:
        raise ValueError(f"Run directories not found: {missing}")

    loaded: list[tuple[str, MetropolisParams, GaugeObservableParams, GradientFlowParams]] = []
    for path in normalized:
        yaml_path = Path(path) / "input.yaml"
        if not yaml_path.exists():
            raise ValueError(f"Missing input.yaml in {path}")
        metro, gauge = load_params(str(yaml_path))
        gradient_flow = load_gradient_flow_params(str(yaml_path))
        loaded.append((path, metro, gauge, gradient_flow))

    ref_path, ref_metro, ref_gauge, ref_gradient_flow = loaded[0]
    ref_metro_dict = asdict(ref_metro)
    ref_gauge_dict = asdict(ref_gauge)
    ref_gradient_flow_dict = asdict(ref_gradient_flow)
    ref_required_metro_cmp = {k: v for k, v in ref_metro_dict.items() if k in REQUIRED_METRO_FIELDS}
    ref_gradient_flow_cmp: dict[str, Any] = {}
    if require_fixed_dt:
        ref_gradient_flow_cmp["dt"] = ref_gradient_flow_dict.get("dt")

    metro_values_by_run = {ref_path: ref_metro_dict}
    gauge_values_by_run = {ref_path: ref_gauge_dict}
    gradient_flow_values_by_run = {ref_path: ref_gradient_flow_dict}

    mismatches: list[str] = []
    for path, metro, gauge, gradient_flow in loaded[1:]:
        current_metro_dict = asdict(metro)
        current_gauge_dict = asdict(gauge)
        current_gradient_flow_dict = asdict(gradient_flow)
        metro_values_by_run[path] = current_metro_dict
        gauge_values_by_run[path] = current_gauge_dict
        gradient_flow_values_by_run[path] = current_gradient_flow_dict
        metro_cmp = {k: v for k, v in current_metro_dict.items() if k in REQUIRED_METRO_FIELDS}
        gradient_flow_cmp: dict[str, Any] = {}
        if require_fixed_dt:
            gradient_flow_cmp["dt"] = current_gradient_flow_dict.get("dt")
        metro_diff = compare_dataclass_dicts(ref_required_metro_cmp, metro_cmp)
        gauge_diff: dict[str, tuple[Any, Any]] = {}
        gradient_flow_diff = compare_dataclass_dicts(ref_gradient_flow_cmp, gradient_flow_cmp)
        if not metro_diff and not gauge_diff and not gradient_flow_diff:
            continue

        parts = [f"Compatibility mismatch for {path} compared to {ref_path}:"]
        for field, (expected, got) in metro_diff.items():
            parts.append(f"  MetropolisParams.{field}: expected {expected!r}, got {got!r}")
        for field, (expected, got) in gauge_diff.items():
            parts.append(f"  GaugeObservableParams.{field}: expected {expected!r}, got {got!r}")
        for field, (expected, got) in gradient_flow_diff.items():
            parts.append(f"  GradientFlowParams.{field}: expected {expected!r}, got {got!r}")
        mismatches.append("\n".join(parts))

    if mismatches:
        raise ValueError("\n".join(mismatches))

    metro_field_names = set(ref_metro_dict)
    gauge_field_names = set(ref_gauge_dict)
    gradient_flow_field_names = set(ref_gradient_flow_dict)
    common_metropolis_fields = _common_field_summary(metro_values_by_run, metro_field_names)
    varying_metropolis_fields = _varying_field_summary(metro_values_by_run, metro_field_names)
    common_gauge_fields = _common_field_summary(gauge_values_by_run, gauge_field_names)
    varying_gauge_fields = _varying_field_summary(gauge_values_by_run, gauge_field_names)
    common_gradient_flow_fields = _common_field_summary(
        gradient_flow_values_by_run,
        gradient_flow_field_names,
    )
    varying_gradient_flow_fields = _varying_field_summary(
        gradient_flow_values_by_run,
        gradient_flow_field_names,
    )
    gradient_flow_dt_by_run = {
        run_dir: values.get("dt") for run_dir, values in gradient_flow_values_by_run.items()
    }
    gradient_flow_t_values_by_run = {
        run_dir: _jsonable_group_value(values.get("t_values"))
        for run_dir, values in gradient_flow_values_by_run.items()
    }
    summary = {
        "reference_run": ref_path,
        "ignored_metropolis_fields": sorted(IGNORED_METRO_FIELDS),
        "required_metropolis_fields": sorted(REQUIRED_METRO_FIELDS),
        "allowed_gradient_flow_differences": sorted(ALLOWED_GRADIENT_FLOW_FIELDS),
        "require_fixed_dt": bool(require_fixed_dt),
        "metropolis_common": common_metropolis_fields,
        "metropolis_varying": varying_metropolis_fields,
        "gauge_common": common_gauge_fields,
        "gauge_varying": varying_gauge_fields,
        "gradient_flow_common": common_gradient_flow_fields,
        "gradient_flow_varying": varying_gradient_flow_fields,
        "gradient_flow_dt_values": _unique_values_for_runs(gradient_flow_dt_by_run),
        "gradient_flow_dt_by_run": gradient_flow_dt_by_run,
        "gradient_flow_t_values_by_run": gradient_flow_t_values_by_run,
    }
    return ref_metro, ref_gauge, summary


class FinalizedAnalysisRunner:
    def __init__(
        self,
        run_dirs: list[str],
        output_root: Path,
        plot_mode: str,
        load_workers: int,
        calc_workers: int | None,
        n_bootstrap: int,
        target_force: float,
        block_size_scan_values: list[int] | None = None,
        block_size_scan_min: int = DEFAULT_BLOCK_SIZE_SCAN_MIN,
        block_size_scan_max: int | None = DEFAULT_BLOCK_SIZE_SCAN_MAX,
        block_size_scan_step: int = DEFAULT_BLOCK_SIZE_SCAN_STEP,
        max_r: int | None = None,
        wilson_flow_time: float | None = None,
        require_fixed_dt: bool = False,
        open_plots: bool = True,
        open_plot_func: Callable[[Path], bool] = open_html_plot,
        input_func: Callable[[str], str] = input,
        print_func: Callable[..., None] = print,
    ):
        if plot_mode != "html":
            raise ValueError("Only --plot-mode html is currently supported.")

        self.run_dirs = [os.path.abspath(path) for path in run_dirs]
        self.output_root = output_root.resolve()
        self.plot_mode = plot_mode
        self.load_workers = int(load_workers)
        self.calc_workers = calc_workers
        self.n_bootstrap = int(n_bootstrap)
        self.target_force = float(target_force)
        self.block_size_scan_min = int(block_size_scan_min)
        self.block_size_scan_max = int(block_size_scan_max) if block_size_scan_max is not None else None
        self.block_size_scan_step = int(block_size_scan_step)
        if self.block_size_scan_min < 1:
            raise ValueError("--min-block-size must be at least 1")
        if self.block_size_scan_step < 1:
            raise ValueError("--block-step must be at least 1")
        if self.block_size_scan_max is not None and self.block_size_scan_max < self.block_size_scan_min:
            raise ValueError("--max-block-size must be greater than or equal to --min-block-size")
        self.max_r = int(max_r) if max_r is not None else None
        if self.max_r is not None and self.max_r < 0:
            raise ValueError("--max-r must be non-negative")
        self.wilson_flow_time = data_organizer._normalize_flow_time(wilson_flow_time)
        self.require_fixed_dt = bool(require_fixed_dt)
        self.block_size_scan_values = (
            sorted({int(value) for value in block_size_scan_values if int(value) >= 1})
            if block_size_scan_values is not None
            else None
        )
        if block_size_scan_values is not None and not self.block_size_scan_values:
            raise ValueError("At least one positive bootstrap block size is required.")
        self.open_plots = bool(open_plots)
        self.open_plot_func = open_plot_func
        self.input = input_func
        self.print = print_func

        self.metro, self.gauge, self.compatibility_summary = validate_run_directories(
            self.run_dirs,
            require_fixed_dt=self.require_fixed_dt,
        )
        self._print_input_compatibility_summary()
        self.selection_preview_dir = (self.output_root / "_selection_previews" / preview_hash(self.run_dirs)).resolve()
        self.selection_preview_dir.mkdir(parents=True, exist_ok=True)
        self.thermalization_selection_path = self.selection_preview_dir / "thermalization_cuts.json"
        self.thermalization_preview_records_by_run: dict[str, list[dict[str, Any]]] = {}
        self.thermalization_preview_path: Path | None = None
        self.thermalization_steps_by_run = self._prompt_thermalization_by_run()
        self.thermalization_steps = common_thermalization_step(self.thermalization_steps_by_run)
        if self.thermalization_steps is not None:
            self.print(f"Thermalization cut set to: {self.thermalization_steps} for all runs")
        else:
            self.print("Thermalization cuts saved per run.")
        self.analysis_id = analysis_hash(
            self.run_dirs,
            self.thermalization_steps_by_run,
            wilson_flow_time=self.wilson_flow_time,
        )
        self.analysis_dir = self._analysis_dir_path()

        self.manifest_path = self.analysis_dir / "manifest.json"
        self.input_runs_path = self.analysis_dir / "input_runs.json"
        self.scan_dir = self.analysis_dir / "scan_cache"
        self.plots_dir = self.analysis_dir / "plots"
        self.wilson_dir = self.analysis_dir / "wilsonloop"
        self.wilson_bootstrap_dir = self.wilson_dir / "bootstrap"
        self.r0_dir = self.analysis_dir / "r0"
        self.r0_bootstrap_dir = self.r0_dir / "bootstrap"
        self.derived_dir = self.analysis_dir / "derived"
        self.derived_bootstrap_dir = self.derived_dir / "bootstrap"
        self.gradient_flow_dir = self.analysis_dir / "gradient_flow"
        self.gradient_flow_bootstrap_dir = self.gradient_flow_dir / "bootstrap"
        self.creutz_dir = self.analysis_dir / "creutz"
        self.creutz_bootstrap_dir = self.creutz_dir / "bootstrap"
        self.bootstrap_block_size_choice_path = self.analysis_dir / "bootstrap_block_size_choice.json"

        self.manifest = self._load_or_init_manifest()
        self._ensure_output_layout()
        self._save_final_thermalization_preview()
        self._write_input_runs()

        self.combined_w_temp = None
        self.combined_gradient_flow = None
        self.gradient_flow_metadata: dict[str, Any] = {}
        self.aggregation = None
        self.calc: Calculator | None = None
        self.unique_rs: list[int] = []
        self.unique_ts: list[int] = []
        self.windows: list[tuple[int, int | None]] = []

    def _analysis_dir_path(self) -> Path:
        flow_part = ""
        if self.wilson_flow_time is not None:
            flow_part = f"__tf_{format_token(self.wilson_flow_time)}"
        folder_name = (
            f"beta_{format_token(self.metro.beta)}"
            f"__L_{self.metro.L0}x{self.metro.L1}x{self.metro.L2}x{self.metro.L3}"
            f"__eps1_{format_token(self.metro.epsilon1)}"
            f"{flow_part}"
            f"__nrun_{len(self.run_dirs)}"
            f"__{self.analysis_id}"
        )
        return self.output_root / folder_name

    def _print_input_compatibility_summary(self) -> None:
        dt_values = self.compatibility_summary.get("gradient_flow_dt_values", [])
        if dt_values:
            dt_text = ", ".join(format_token(value) for value in dt_values)
            mode = "fixed" if len(dt_values) == 1 else "varies"
            self.print(f"Gradient-flow dt ({mode}): {dt_text}")
        t_values_by_run = self.compatibility_summary.get("gradient_flow_t_values_by_run", {})
        unique_t_value_sets = _unique_values_for_runs(t_values_by_run) if t_values_by_run else []
        if unique_t_value_sets:
            self.print(
                "Gradient-flow t_values set(s): "
                f"{len(unique_t_value_sets)} unique list(s); combined summaries use all observed flow times."
            )
        varying_metro = self.compatibility_summary.get("metropolis_varying", {})
        if "nSweep" in varying_metro:
            values = varying_metro["nSweep"].get("unique_values", [])
            self.print(f"nSweep varies across runs: {', '.join(format_token(value) for value in values)}")

    def _load_saved_thermalization_cuts(self) -> dict[str, int]:
        payload = load_json(self.thermalization_selection_path, default={})
        raw_cuts = payload.get("cuts_by_run", payload) if isinstance(payload, dict) else {}
        if not isinstance(raw_cuts, dict):
            return {}

        cuts_by_run: dict[str, int] = {}
        for run_dir in self.run_dirs:
            value = raw_cuts.get(run_dir)
            if value is None:
                value = raw_cuts.get(os.path.abspath(run_dir))
            if value is None:
                continue
            try:
                cut = int(value)
            except (TypeError, ValueError):
                continue
            if cut < 0:
                continue
            cuts_by_run[os.path.abspath(run_dir)] = cut
        return cuts_by_run

    def _save_thermalization_selection_state(self, cuts_by_run: dict[str, int]) -> None:
        normalized = normalize_thermalization_steps_by_run(
            [run_dir for run_dir in self.run_dirs if run_dir in cuts_by_run],
            cuts_by_run,
        )
        save_json(
            self.thermalization_selection_path,
            {
                "schema_version": SCHEMA_VERSION,
                "run_dirs": self.run_dirs,
                "cuts_by_run": normalized,
                "updated_at": utc_now(),
            },
        )

    def _thermalization_preview_records_for_run(self, run_dir: str) -> list[dict[str, Any]]:
        run_dir = os.path.abspath(run_dir)
        cached = self.thermalization_preview_records_by_run.get(run_dir)
        if cached is None:
            cached = build_thermalization_preview([run_dir], include_combined=False)
            self.thermalization_preview_records_by_run[run_dir] = cached
        return cached

    def _selection_preview_path_for_run(self, run_dir: str) -> Path:
        return self.selection_preview_dir / thermalization_preview_filename(run_dir)

    def _analysis_preview_path_for_run(self, run_dir: str) -> Path:
        return self.plots_dir / "thermalization_by_run" / thermalization_preview_filename(run_dir)

    def _save_selection_thermalization_preview(self, run_dir: str, suggested_cut: int) -> Path | None:
        records = self._thermalization_preview_records_for_run(run_dir)
        if not records:
            return None

        preview_path = self._selection_preview_path_for_run(run_dir)
        save_thermalization_plot(
            preview_path,
            records,
            suggested_cut=int(suggested_cut),
        )
        self.thermalization_preview_path = preview_path
        return preview_path

    def _open_plot(self, path: Path) -> bool:
        if not self.open_plots:
            return False
        try:
            return bool(self.open_plot_func(path))
        except Exception:
            return False

    def _prompt_single_run_thermalization(self, run_dir: str, suggested_cut: int) -> tuple[int, bool]:
        run_name = Path(run_dir).name
        while True:
            response = self.input(
                f"Thermalization cut for {run_name} is {int(suggested_cut)}. "
                "Press Enter to accept, type a new integer, or append * to apply to all runs: "
            ).strip()
            if response == "" or response.lower() in {"y", "yes"}:
                return int(suggested_cut), False
            apply_to_all = response.endswith("*")
            if apply_to_all:
                response = response[:-1].strip()
            try:
                value = int(response)
            except ValueError:
                self.print(f"Invalid thermalization cut: {response!r}")
                continue
            if value < 0:
                self.print("Thermalization cut must be non-negative.")
                continue
            return value, apply_to_all

    def _prompt_thermalization_by_run(self) -> dict[str, int]:
        default_cut = int(run_evaluation.THERMALIZATION_STEPS)
        self.print("Selected run directories:")
        for index, path in enumerate(self.run_dirs, start=1):
            self.print(f"  {index}. {path}")

        cuts_by_run = self._load_saved_thermalization_cuts()
        if cuts_by_run:
            self.print(f"Loaded saved thermalization cuts for {len(cuts_by_run)} run(s).")

        for index, run_dir in enumerate(self.run_dirs, start=1):
            if run_dir in cuts_by_run:
                self.print(
                    f"Reusing thermalization cut {cuts_by_run[run_dir]} for run {index}/{len(self.run_dirs)}: {run_dir}"
                )
                continue

            self.print(f"\nRun {index}/{len(self.run_dirs)}: {run_dir}")
            preview_path = self._save_selection_thermalization_preview(run_dir, default_cut)
            if preview_path is not None:
                opened = self._open_plot(preview_path)
                action = "Opened" if opened else "Saved"
                self.print(f"{action} thermalization preview: {preview_path}")
            else:
                self.print("No thermalization preview data was available for this run.")

            cut, apply_to_all = self._prompt_single_run_thermalization(run_dir, default_cut)
            if apply_to_all:
                cuts_by_run = {path: int(cut) for path in self.run_dirs}
                self._save_thermalization_selection_state(cuts_by_run)
                self.print(f"Saved thermalization cut {cut} for all runs.")
                break
            cuts_by_run[run_dir] = int(cut)
            self._save_thermalization_selection_state(cuts_by_run)
            self.print(f"Saved thermalization cut {cut} for {Path(run_dir).name}.")

        normalized = normalize_thermalization_steps_by_run(self.run_dirs, cuts_by_run)
        self.print("Thermalization summary:")
        for run_dir, cut in normalized.items():
            self.print(f"  {Path(run_dir).name}: {cut}")
        return normalized

    def _write_thermalization_index(self, path: Path, entries: list[dict[str, Any]]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        rows: list[str] = []
        for entry in entries:
            run_dir = str(entry["run_dir"])
            run_name = html.escape(Path(run_dir).name or run_dir)
            cut = int(entry["cut"])
            plot_path = entry.get("plot_path")
            link_html = "not available"
            if plot_path is not None:
                rel_path = os.path.relpath(str(plot_path), start=str(path.parent))
                link_html = f'<a href="{html.escape(rel_path)}">interactive preview</a>'
            rows.append(
                "<tr>"
                f"<td><code>{run_name}</code></td>"
                f"<td><code>{html.escape(run_dir)}</code></td>"
                f"<td>{cut}</td>"
                f"<td>{link_html}</td>"
                "</tr>"
            )

        body_rows = "\n".join(rows) if rows else "<tr><td colspan='4'>No thermalization previews available.</td></tr>"
        document = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Thermalization Selection Summary</title>
  <style>
    body {{ font-family: sans-serif; margin: 2rem; color: #222; }}
    table {{ border-collapse: collapse; width: 100%; }}
    th, td {{ border: 1px solid #ddd; padding: 0.6rem; text-align: left; vertical-align: top; }}
    th {{ background: #f5f5f5; }}
    code {{ font-family: ui-monospace, SFMono-Regular, Menlo, monospace; }}
  </style>
</head>
<body>
  <h1>Thermalization Selection Summary</h1>
  <p>Each run was reviewed individually before the combined analysis was built.</p>
  <table>
    <thead>
      <tr>
        <th>Run</th>
        <th>Directory</th>
        <th>Cut</th>
        <th>Preview</th>
      </tr>
    </thead>
    <tbody>
      {body_rows}
    </tbody>
  </table>
</body>
</html>
"""
        path.write_text(document, encoding="utf-8")

    def _save_final_thermalization_preview(self) -> None:
        entries: list[dict[str, Any]] = []
        for run_dir in self.run_dirs:
            records = self._thermalization_preview_records_for_run(run_dir)
            plot_path: Path | None = None
            if records:
                plot_path = self._analysis_preview_path_for_run(run_dir)
                save_thermalization_plot(
                    plot_path,
                    records,
                    suggested_cut=int(self.thermalization_steps_by_run[run_dir]),
                )
            entries.append(
                {
                    "run_dir": run_dir,
                    "cut": int(self.thermalization_steps_by_run[run_dir]),
                    "plot_path": plot_path,
                }
            )

        final_path = self.plots_dir / "thermalization_preview.html"
        self._write_thermalization_index(final_path, entries)
        self.manifest["thermalization_preview_plot"] = str(final_path)
        self.manifest["thermalization_preview_plots_by_run"] = {
            str(entry["run_dir"]): str(entry["plot_path"])
            for entry in entries
            if entry["plot_path"] is not None
        }
        self._save_manifest()

    def _load_or_init_manifest(self) -> dict[str, Any]:
        manifest = load_json(self.manifest_path, default=None)
        if manifest is None:
            return {
                "schema_version": SCHEMA_VERSION,
                "analysis_id": self.analysis_id,
                "created_at": utc_now(),
                "updated_at": utc_now(),
                "status": {
                    "scan_cache_built": False,
                    "wilsonloop_complete": False,
                    "r0_complete": False,
                    "derived_complete": False,
                    "gradient_flow_complete": False,
                    "creutz_complete": False,
                },
                "run_evaluation_version": run_evaluation.CALC_VERSION,
                "calculator_helper": "fit_r0_from_potential_data",
                "thermalization_mode": "per_run",
                "thermalization_steps": self.thermalization_steps,
                "thermalization_steps_by_run": self.thermalization_steps_by_run,
                "wilson_flow_time": None if self.wilson_flow_time is None else float(self.wilson_flow_time),
                "n_bootstrap": int(self.n_bootstrap),
                "target_force": float(self.target_force),
                "block_size": None,
                "block_size_scan_values": self.block_size_scan_values,
                "block_size_scan_min": self.block_size_scan_min,
                "block_size_scan_max": self.block_size_scan_max,
                "block_size_scan_step": self.block_size_scan_step,
                "plot_mode": self.plot_mode,
                "calc_workers": self.calc_workers,
                "load_workers": self.load_workers,
                "compatibility_summary": self.compatibility_summary,
            }

        stored_cuts_raw = manifest.get("thermalization_steps_by_run")
        if stored_cuts_raw is None:
            stored_step = manifest.get("thermalization_steps")
            if stored_step is None:
                raise ValueError("Existing analysis is missing thermalization cut information.")
            stored_cuts = {run_dir: int(stored_step) for run_dir in self.run_dirs}
        else:
            stored_cuts = normalize_thermalization_steps_by_run(self.run_dirs, stored_cuts_raw)

        if stored_cuts != self.thermalization_steps_by_run:
            raise ValueError(
                "Existing analysis uses different per-run thermalization cuts than this run requested."
            )
        if int(manifest.get("n_bootstrap")) != int(self.n_bootstrap):
            raise ValueError(
                f"Existing analysis uses n_bootstrap={manifest.get('n_bootstrap')} "
                f"but this run requested {self.n_bootstrap}."
            )
        if not np.isclose(float(manifest.get("target_force")), float(self.target_force)):
            raise ValueError(
                f"Existing analysis uses target_force={manifest.get('target_force')} "
                f"but this run requested {self.target_force}."
            )
        stored_flow = data_organizer._normalize_flow_time(manifest.get("wilson_flow_time"))
        if stored_flow != self.wilson_flow_time:
            raise ValueError(
                f"Existing analysis uses wilson_flow_time={manifest.get('wilson_flow_time')} "
                f"but this run requested {self.wilson_flow_time}."
            )
        status = manifest.setdefault("status", {})
        status.setdefault("gradient_flow_complete", False)
        status.setdefault("creutz_complete", False)
        manifest["updated_at"] = utc_now()
        return manifest

    def _ensure_output_layout(self) -> None:
        for path in [
            self.analysis_dir,
            self.scan_dir,
            self.plots_dir,
            self.wilson_dir,
            self.wilson_bootstrap_dir,
            self.r0_dir,
            self.r0_bootstrap_dir,
            self.derived_dir,
            self.derived_bootstrap_dir,
            self.gradient_flow_dir,
            self.gradient_flow_bootstrap_dir,
            self.creutz_dir,
            self.creutz_bootstrap_dir,
        ]:
            path.mkdir(parents=True, exist_ok=True)
        self._save_manifest()

    def _write_input_runs(self) -> None:
        payload = {
            "run_dirs": self.run_dirs,
            "thermalization_steps_by_run": self.thermalization_steps_by_run,
            "wilson_flow_time": None if self.wilson_flow_time is None else float(self.wilson_flow_time),
            "compatibility_summary": self.compatibility_summary,
        }
        save_json(self.input_runs_path, payload)

    def _save_manifest(self) -> None:
        self.manifest["updated_at"] = utc_now()
        save_json(self.manifest_path, self.manifest)

    def _saved_bootstrap_block_size(self) -> int | None:
        choice_payload = load_json(self.bootstrap_block_size_choice_path, default=None)
        if isinstance(choice_payload, dict):
            try:
                value = int(choice_payload.get("block_size"))
            except (TypeError, ValueError):
                value = 0
            if value >= 1:
                return value

        try:
            value = int(self.manifest.get("block_size"))
        except (TypeError, ValueError):
            value = 0
        return value if value >= 1 else None

    def _compute_tau_hint(self, file_data) -> float | None:
        if not any(obs.name in {"plaquette", "retrace"} for obs in getattr(file_data, "observables", [])):
            return None

        calc = Calculator(file_data, n_bootstrap=self.n_bootstrap, step_size=1)
        for obs_name in ("plaquette", "retrace"):
            try:
                tau = calc.get_variable("tau_int", obs_name=obs_name).get()
            except KeyError:
                continue
            if tau is not None and np.isfinite(tau):
                return float(tau)
        return None

    def _bootstrap_scan_total_length(self, file_data) -> int:
        n_configurations = getattr(file_data, "n_configurations", None)
        try:
            if n_configurations is not None and int(n_configurations) > 0:
                return int(n_configurations)
        except (TypeError, ValueError):
            pass

        wilson_by_pair = getattr(file_data, "wilson_by_pair", None)
        if wilson_by_pair:
            first_series = next(iter(wilson_by_pair.values()), None)
            if first_series is not None:
                return int(len(first_series))
        wilson_by_flow_pair = getattr(file_data, "wilson_by_flow_pair", None)
        if wilson_by_flow_pair:
            selected_flow = self.wilson_flow_time
            first_series = next(
                (
                    series
                    for (flow_time, _r_val, _t_val), series in wilson_by_flow_pair.items()
                    if flow_time == selected_flow
                ),
                None,
            )
            if first_series is not None:
                return int(len(first_series))

        lengths = [len(obs.values) for obs in getattr(file_data, "observables", [])]
        return int(min(lengths)) if lengths else 0

    def _resolved_bootstrap_block_sizes(self, file_data) -> list[int]:
        if self.block_size_scan_values is not None:
            return list(self.block_size_scan_values)

        total_length = self._bootstrap_scan_total_length(file_data)
        if total_length <= 0:
            raise ValueError("Cannot infer bootstrap block-size scan range from empty data.")

        max_block = self.block_size_scan_max
        if max_block is None:
            max_block = max(
                self.block_size_scan_min,
                int(ceil(total_length / DEFAULT_BLOCK_SIZE_SCAN_MAX_FRACTION)),
            )

        values = list(range(self.block_size_scan_min, int(max_block) + 1, self.block_size_scan_step))
        if values and values[-1] != int(max_block):
            values.append(int(max_block))

        block_sizes = sorted({int(value) for value in values if int(value) >= 1})
        if not block_sizes:
            raise ValueError("No valid bootstrap block sizes selected.")
        return block_sizes

    def _prompt_bootstrap_block_size(self, recommended_block_size: int) -> int:
        while True:
            response = self.input(
                "Bootstrap block size in saved configurations "
                f"is {int(recommended_block_size)}. Press Enter to accept, or type a new integer: "
            ).strip()
            if response == "" or response.lower() in {"y", "yes"}:
                return int(recommended_block_size)
            try:
                value = int(response)
            except ValueError:
                self.print(f"Invalid bootstrap block size: {response!r}")
                continue
            if value < 1:
                self.print("Bootstrap block size must be at least 1.")
                continue
            return value

    def _select_bootstrap_block_size(self, file_data) -> int:
        saved_block_size = self._saved_bootstrap_block_size()
        if saved_block_size is not None:
            self.print(f"Reusing bootstrap block size: {saved_block_size}")
            self.manifest["block_size"] = int(saved_block_size)
            self._save_manifest()
            return int(saved_block_size)

        self.print("\nBuilding bootstrap block-size diagnostic for representative W(R,T) points...")
        block_size_scan_values = self._resolved_bootstrap_block_sizes(file_data)
        self.manifest["block_size_scan_values"] = block_size_scan_values
        self._save_manifest()
        scan_path = self.scan_dir / "bootstrap_block_size_scan.json"
        plot_path = self.plots_dir / "bootstrap_block_size.html"
        scan_records, selected_pairs = build_bootstrap_block_size_scan(
            file_data,
            block_size_scan_values,
            n_bootstrap=self.n_bootstrap,
            flow_time=self.wilson_flow_time,
        )
        tau_hint = self._compute_tau_hint(file_data)
        recommended = recommend_bootstrap_block_size(
            scan_records,
            fallback=run_evaluation.DEFAULT_BLOCK_SIZE,
        )

        save_json(
            scan_path,
            {
                "schema_version": SCHEMA_VERSION,
                "block_size_scan_values": block_size_scan_values,
                "block_size_scan_total_length": self._bootstrap_scan_total_length(file_data),
                "wilson_flow_time": None if self.wilson_flow_time is None else float(self.wilson_flow_time),
                "selected_pairs": [{"R": int(r_val), "T": int(t_val)} for r_val, t_val in selected_pairs],
                "tau_int": tau_hint,
                "recommended_block_size": int(recommended),
                "records": scan_records,
                "saved_at": utc_now(),
            },
        )
        save_bootstrap_block_size_plot(
            plot_path,
            scan_records,
            recommended_block_size=int(recommended),
        )
        opened = self._open_plot(plot_path)
        action = "Opened" if opened else "Saved"
        self.print(f"{action} bootstrap block-size plot: {plot_path}")
        if tau_hint is None:
            self.print("tau_int unavailable from loaded data; using the error-vs-block-size diagnostic.")
        else:
            self.print(f"tau_int hint: {tau_hint:.6g}")

        block_size = self._prompt_bootstrap_block_size(int(recommended))
        save_json(
            self.bootstrap_block_size_choice_path,
            {
                "schema_version": SCHEMA_VERSION,
                "block_size": int(block_size),
                "recommended_block_size": int(recommended),
                "scan_path": str(scan_path),
                "plot_path": str(plot_path),
                "wilson_flow_time": None if self.wilson_flow_time is None else float(self.wilson_flow_time),
                "tau_int": tau_hint,
                "saved_at": utc_now(),
            },
        )
        self.manifest["block_size"] = int(block_size)
        self.manifest["bootstrap_block_size_scan"] = str(scan_path)
        self.manifest["bootstrap_block_size_plot"] = str(plot_path)
        self.manifest["tau_int"] = tau_hint
        self._save_manifest()
        return int(block_size)

    def _load_data(self) -> None:
        if self.combined_w_temp is not None and self.calc is not None:
            return

        combined_w_temp, aggregation = run_evaluation._load_combined_w_temp_filtered(
            self.run_dirs,
            pair_filter=wilson_flow_time_pair_filter(self.wilson_flow_time),
            filter_label=wilson_flow_time_filter_label(self.wilson_flow_time),
            verbose=True,
            prefix="[finalize_analysis]",
            load_workers=self.load_workers,
            thermalization_steps_by_run=self.thermalization_steps_by_run,
        )
        combined_gradient_flow, gradient_flow_metadata = run_evaluation._load_combined_gradient_flow_obs(
            self.run_dirs,
            load_workers=self.load_workers,
            thermalization_steps_by_run=self.thermalization_steps_by_run,
        )
        if combined_w_temp is None:
            raise RuntimeError("No combined W_temp data could be loaded from the selected runs.")
        
        self.print(f"Combined W_temp data loaded with aggregation: {aggregation}")

        block_size = self._select_bootstrap_block_size(combined_w_temp)
        self.combined_w_temp = combined_w_temp
        self.combined_gradient_flow = combined_gradient_flow
        self.gradient_flow_metadata = gradient_flow_metadata
        self.aggregation = aggregation
        self.calc = Calculator(
            combined_w_temp,
            n_bootstrap=self.n_bootstrap,
            step_size=block_size,
        )
        available_flow_times = self.calc.get_available_flow_times()
        if self.wilson_flow_time not in available_flow_times:
            available_text = ", ".join(
                "unflowed" if value is None else f"{float(value):g}"
                for value in available_flow_times
            ) or "none"
            requested_text = "unflowed" if self.wilson_flow_time is None else f"{float(self.wilson_flow_time):g}"
            raise RuntimeError(
                f"Requested Wilson-loop flow time {requested_text} is not available. "
                f"Available Wilson-loop flow times: {available_text}."
            )
        if self.wilson_flow_time is not None:
            self.print(f"Using Wilson-loop flow time t_f={self.wilson_flow_time:g}.")
        self.unique_rs = [int(r) for r in self.calc.get_unique_Rs(flow_time=self.wilson_flow_time)]
        self.unique_ts = [int(t) for t in self.calc.get_unique_Ts(flow_time=self.wilson_flow_time)]
        if not self.unique_rs or not self.unique_ts:
            requested_text = "unflowed" if self.wilson_flow_time is None else f"{float(self.wilson_flow_time):g}"
            raise RuntimeError(f"No Wilson-loop R/T pairs are available for flow time {requested_text}.")
        self.windows = enumerate_t_windows(self.unique_ts)

        self.manifest["aggregation"] = aggregation
        self.manifest["block_size"] = block_size
        self.manifest["available_R"] = self.unique_rs
        self.manifest["available_T"] = self.unique_ts
        self.manifest["wilson_flow_time"] = None if self.wilson_flow_time is None else float(self.wilson_flow_time)
        self.manifest["available_wilson_flow_times"] = [
            None if value is None else float(value)
            for value in available_flow_times
        ]
        self.manifest["gradient_flow_metadata"] = gradient_flow_metadata
        self._save_manifest()

    def _build_scan_cache(self) -> None:
        self._load_data()
        self.print("\nBuilding scan cache for V(R) analysis...")
        assert self.calc is not None

        wrt_path = self.scan_dir / "wrt_scan.json"
        eff_path = self.scan_dir / "effective_mass_scan.json"
        vr_path = self.scan_dir / "v_r_scan_candidates.json"

        if not wrt_path.exists():
            save_json(wrt_path, build_wrt_scan(self.calc, self.unique_rs, self.unique_ts, flow_time=self.wilson_flow_time))
        if not eff_path.exists():
            save_json(
                eff_path,
                build_effective_mass_scan(
                    self.calc,
                    self.unique_rs,
                    self.unique_ts,
                    flow_time=self.wilson_flow_time,
                ),
            )
        if not vr_path.exists():
            save_json(vr_path, build_v_r_scan(self.calc, self.unique_rs, self.windows, flow_time=self.wilson_flow_time))

        self.manifest["status"]["scan_cache_built"] = True
        self._save_manifest()

    def _load_v_choices(self) -> dict[str, Any]:
        return load_json(self.wilson_dir / "choices.json", default={"choices": {}})

    def _load_v_summary(self) -> dict[str, Any]:
        return load_json(self.wilson_dir / "V_R_summary.json", default={"results": {}})

    def _save_v_state(self, choices_payload: dict[str, Any], summary_payload: dict[str, Any]) -> None:
        save_json(self.wilson_dir / "choices.json", choices_payload)
        save_json(self.wilson_dir / "V_R_summary.json", summary_payload)

    def _ignored_r_ge(self, choices_payload: dict[str, Any]) -> int | None:
        ignored = choices_payload.get("ignored_R", {})
        if isinstance(ignored, dict) and ignored.get("R_ge") is not None:
            return int(ignored["R_ge"])
        if choices_payload.get("ignore_R_ge") is not None:
            return int(choices_payload["ignore_R_ge"])
        return None

    def _is_v_r_ignored(self, r_value: int, choices_payload: dict[str, Any]) -> bool:
        if self.max_r is not None and int(r_value) > self.max_r:
            return True
        ignored_r_ge = self._ignored_r_ge(choices_payload)
        return ignored_r_ge is not None and int(r_value) >= ignored_r_ge

    def _save_v_stop(
        self,
        choices_payload: dict[str, Any],
        summary_payload: dict[str, Any],
        stop_before_r: int,
        *,
        source: str,
    ) -> None:
        kept_rs = [
            int(key)
            for key in choices_payload.get("choices", {})
            if int(key) < int(stop_before_r)
        ]
        max_r = max(kept_rs) if kept_rs else int(stop_before_r) - 1
        choices_payload["ignored_R"] = {
            "R_ge": int(stop_before_r),
            "max_R": int(max_r),
            "source": source,
            "saved_at": utc_now(),
        }
        # Keep the flat fields for readability and compatibility with older output snapshots.
        choices_payload["ignore_R_ge"] = int(stop_before_r)
        choices_payload["max_R"] = int(max_r)
        self._save_v_state(choices_payload, summary_payload)

    def _candidate_windows_for_r(self, r_value: int) -> list[dict[str, Any]]:
        v_scan = load_json(self.scan_dir / "v_r_scan_candidates.json", default=[])
        return [
            row
            for row in v_scan
            if int(row["R"]) == int(r_value)
        ]

    def _prompt_window_choice(
        self,
        r_value: int,
        candidates: list[dict[str, Any]],
        previous_choice: tuple[int, int | None] | None = None,
    ) -> tuple[int, int | None] | None:
        if not candidates:
            self.print(f"\nR = {r_value}")
            self.print(f"No finite V(R) candidate windows are available for R={r_value}.")
            while True:
                response = self.input(
                    f"Type `stop` to ignore R >= {r_value} and continue, or interrupt to abort: "
                ).strip().lower()
                if response == "stop":
                    return None
                self.print("No fit window can be selected for this R.")

        self.print(f"\nR = {r_value}")
        while True:
            prompt = "Enter t_min and optional t_max as `t_min` or `t_min,t_max` (use `none` for no t_max)"
            if previous_choice is not None:
                prompt += f", or press Enter to reuse {window_label(*previous_choice)}"
            prompt += f"; type `stop` to ignore R >= {r_value} and continue"
            prompt += ": "
            response = self.input(prompt).strip()
            if response.lower() == "stop":
                return None
            if not response:
                if previous_choice is None:
                    self.print("A t_min selection is required for the first R.")
                    continue
                t_min, t_max = previous_choice
            else:
                try:
                    parts = [part.strip().lower() for part in response.split(",")]
                    if len(parts) == 1:
                        t_min = int(parts[0])
                        t_max = None
                    elif len(parts) == 2:
                        t_min = int(parts[0])
                        t_max = None if parts[1] in {"", "none"} else int(parts[1])
                    else:
                        raise ValueError
                except ValueError:
                    self.print("Please enter either `t_min` or `t_min,t_max`.")
                    continue

            for row in candidates:
                if int(row["t_min"]) == t_min and row["t_max"] == t_max:
                    return t_min, t_max

            self.print(f"Window {window_label(t_min, t_max)} is not available for R={r_value}.")

    def _finalize_wilsonloop(self) -> None:
        self._build_scan_cache()
        self.print("\nScan cache built. Proceeding to finalize V(R) results for each R...")
        assert self.calc is not None

        wrt_records = load_json(self.scan_dir / "wrt_scan.json", default=[])
        effective_mass_records = load_json(self.scan_dir / "effective_mass_scan.json", default=[])
        v_choices = self._load_v_choices()
        v_summary = self._load_v_summary()
        plot_path = self.plots_dir / "effective_mass_all_R.html"
        save_effective_mass_plot(plot_path, wrt_records, effective_mass_records)
        self.print(f"\nSaved effective-mass plot with R dropdown: {plot_path}")
        if self.max_r is not None:
            self.print(f"Limiting V(R) finalization to R <= {self.max_r}.")
        previous_choice: tuple[int, int | None] | None = None

        for r_value in sorted(self.unique_rs):
            key = str(int(r_value))
            if self._is_v_r_ignored(int(r_value), v_choices):
                self._save_v_stop(v_choices, v_summary, int(r_value), source="max_r" if self.max_r is not None else "interactive")
                self.print(f"\nIgnoring R >= {int(r_value)} for V(R).")
                break

            bootstrap_path = self.wilson_bootstrap_dir / f"V_R_R_{key}.npy"
            if key in v_choices["choices"] and bootstrap_path.exists():
                choice = v_choices["choices"][key]
                previous_choice = (int(choice["t_min"]), choice["t_max"])
                continue

            candidates = self._candidate_windows_for_r(r_value)
            window_choice = self._prompt_window_choice(r_value, candidates, previous_choice=previous_choice)
            if window_choice is None:
                self._save_v_stop(v_choices, v_summary, int(r_value), source="interactive")
                self.print(f"Ignoring R >= {int(r_value)} for V(R).")
                break
            t_min, t_max = window_choice
            previous_choice = (t_min, t_max)

            var = self.calc.get_variable(
                "V_R",
                R=int(r_value),
                t_min=int(t_min),
                t_max=t_max,
                flow_time=self.wilson_flow_time,
            )
            value = var.get()
            if value is None or not np.isfinite(value):
                raise RuntimeError(f"Selected window produced no finite V(R) result for R={r_value}.")

            bootstrap_samples = np.asarray(var.bootstrap(), dtype=float)
            np.save(bootstrap_path, bootstrap_samples)

            choice_payload = {
                "R": int(r_value),
                "t_min": int(t_min),
                "t_max": t_max,
                "flow_time": None if self.wilson_flow_time is None else float(self.wilson_flow_time),
                "plot_path": str(plot_path),
                "saved_at": utc_now(),
            }
            summary_payload = {
                "value": float(value),
                "err": float(var.err()) if var.err() is not None else None,
                "flow_time": None if self.wilson_flow_time is None else float(self.wilson_flow_time),
                "fit_C": (
                    float(var.parameters.get("fit_C"))
                    if var.parameters.get("fit_C") is not None and np.isfinite(var.parameters.get("fit_C"))
                    else None
                ),
                "bootstrap_path": str(bootstrap_path),
            }
            v_choices["choices"][key] = choice_payload
            v_summary["results"][key] = summary_payload
            self._save_v_state(v_choices, v_summary)

        self.manifest["status"]["wilsonloop_complete"] = True
        self._save_manifest()

    def _load_selected_v_results(self) -> dict[int, dict[str, Any]]:
        choices_payload = self._load_v_choices()
        summary_payload = self._load_v_summary()
        results: dict[int, dict[str, Any]] = {}
        for key, choice in choices_payload.get("choices", {}).items():
            r_value = int(key)
            if self._is_v_r_ignored(r_value, choices_payload):
                continue
            summary = summary_payload.get("results", {}).get(key)
            if summary is None:
                continue
            bootstrap_path = Path(summary["bootstrap_path"])
            if not bootstrap_path.is_absolute():
                bootstrap_path = (self.analysis_dir / bootstrap_path).resolve()
            if not bootstrap_path.exists():
                bootstrap_path = self.wilson_bootstrap_dir / f"V_R_R_{key}.npy"
            if not bootstrap_path.exists():
                continue
            results[r_value] = {
                "R": r_value,
                "t_min": int(choice["t_min"]),
                "t_max": choice["t_max"],
                "value": float(summary["value"]),
                "err": float(summary["err"]) if summary.get("err") is not None else None,
                "fit_C": summary.get("fit_C"),
                "bootstrap_samples": np.load(bootstrap_path),
                "plot_path": choice["plot_path"],
            }
        return results

    def _prompt_r_min_choice(self, r0_records: list[dict[str, Any]]) -> int:
        if not r0_records:
            raise RuntimeError("No candidate r_min values available for the locked V(R) results.")

        self.print("\nCandidate r_min values:")
        for row in sorted(r0_records, key=lambda item: item["r_min"]):
            chi2_text = ""
            if row.get("chi2_dof") is not None and np.isfinite(row["chi2_dof"]):
                chi2_text = f", chi2/dof={row['chi2_dof']:.3f}"
            self.print(
                f"  r_min={row['r_min']}: r0={row['r0']:.6g}"
                + (f" +/- {row['err']:.3g}" if row.get("err") is not None else "")
                + chi2_text
            )
        valid = {int(row["r_min"]) for row in r0_records}
        while True:
            response = self.input("Choose r_min: ").strip()
            if not response:
                self.print("An r_min selection is required.")
                continue
            try:
                r_min = int(response)
            except ValueError:
                self.print("Please enter an integer r_min.")
                continue
            if r_min not in valid:
                self.print(f"r_min={r_min} is not one of the available candidates: {sorted(valid)}")
                continue
            return r_min

    def _finalize_r0(self) -> None:
        choice_path = self.r0_dir / "choice.json"
        result_path = self.r0_dir / "r0_result.json"
        bootstrap_path = self.r0_bootstrap_dir / "r0.npy"
        selected_v = self._load_selected_v_results()
        if choice_path.exists() and result_path.exists() and bootstrap_path.exists():
            existing_result = load_json(result_path, default={})
            existing_fit_rs = [int(r) for r in existing_result.get("fit_rs", [])]
            if existing_fit_rs and all(r in selected_v for r in existing_fit_rs):
                self.manifest["status"]["r0_complete"] = True
                self._save_manifest()
                return
            self.print("\nExisting r0 result uses V(R) points outside the active cutoff; recomputing.")

        if len(selected_v) < 3:
            raise RuntimeError("At least three finalized V(R) points are required before fitting r0.")

        r0_records, curve_records = build_locked_r0_scan(
            selected_v,
            target_force=self.target_force,
            n_bootstrap=self.n_bootstrap,
        )
        save_json(
            self.scan_dir / "r0_scan_candidates.json",
            {
                "r0_scan": r0_records,
                "cornell_curves": curve_records,
            },
        )

        stability_path = self.plots_dir / "r0_stability.html"
        cornell_path = self.plots_dir / "cornell_by_rmin.html"
        save_r0_stability_plot(stability_path, r0_records)
        save_cornell_plot(cornell_path, curve_records, r0_records)
        self.print(f"\nSaved r0 stability plot: {stability_path}")
        self.print(f"Saved Cornell fit plot: {cornell_path}")

        r_min = self._prompt_r_min_choice(r0_records)
        fit_rs = [r for r in sorted(selected_v) if r >= int(r_min)]
        rs = np.asarray(fit_rs, dtype=float)
        vs = np.asarray([selected_v[r]["value"] for r in fit_rs], dtype=float)
        errs = np.asarray(
            [selected_v[r]["err"] if selected_v[r]["err"] is not None else 1.0 for r in fit_rs],
            dtype=float,
        )
        boot_matrix = np.stack(
            [np.asarray(selected_v[r]["bootstrap_samples"], dtype=float) for r in fit_rs],
            axis=0,
        )
        fit_result = fit_r0_from_potential_data(
            rs,
            vs,
            errs=errs,
            bootstrap_matrix=boot_matrix,
            target_force=self.target_force,
            n_bootstrap=self.n_bootstrap,
        )

        np.save(bootstrap_path, np.asarray(fit_result["bootstrap_samples"], dtype=float))
        save_json(
            choice_path,
            {
                "r_min": int(r_min),
                "flow_time": None if self.wilson_flow_time is None else float(self.wilson_flow_time),
                "plot_paths": {
                    "stability": str(stability_path),
                    "cornell": str(cornell_path),
                },
                "saved_at": utc_now(),
            },
        )
        save_json(
            self.r0_dir / "cornell_points.json",
            [
                selected_v[r] | {"R": int(r)}
                for r in fit_rs
            ],
        )
        save_json(
            result_path,
            {
                "r0": float(fit_result["r0"]),
                "r0_err": float(fit_result["r0_err"]) if fit_result["r0_err"] is not None else None,
                "target_force": float(self.target_force),
                "flow_time": None if self.wilson_flow_time is None else float(self.wilson_flow_time),
                "cornell_params": fit_result["cornell_params"],
                "chi2": float(fit_result["chi2"]) if fit_result["chi2"] is not None else None,
                "dof": int(fit_result["dof"]),
                "chi2_dof": float(fit_result["chi2_dof"]) if fit_result["chi2_dof"] is not None else None,
                "fit_rs": [int(r) for r in fit_rs],
                "bootstrap_path": str(bootstrap_path),
            },
        )

        self.manifest["status"]["r0_complete"] = True
        self._save_manifest()

    def _finalize_derived(self) -> None:
        summary_path = self.derived_dir / "summary.json"
        r0_result = load_json(self.r0_dir / "r0_result.json", default=None)
        if r0_result is None:
            raise RuntimeError("r0 result is missing; cannot compute derived quantities.")
        if summary_path.exists():
            existing_summary = load_json(summary_path, default=None)
            if (
                existing_summary is not None
                and existing_summary.get("r0") is not None
                and np.isclose(float(existing_summary["r0"]), float(r0_result["r0"]))
            ):
                self.manifest["status"]["derived_complete"] = True
                self._save_manifest()
                return
            self.print("\nExisting derived summary is stale relative to the active r0 result; recomputing.")

        r0_bootstrap = np.load(self.r0_bootstrap_dir / "r0.npy")

        r0_val = float(r0_result["r0"])
        r0_err = r0_result.get("r0_err")
        if not np.isfinite(r0_val) or r0_val <= 0:
            raise RuntimeError("Locked r0 is not finite and positive; cannot compute derived quantities.")

        a_boot = np.full_like(r0_bootstrap, np.nan, dtype=float)
        valid_r0 = np.isfinite(r0_bootstrap) & (r0_bootstrap > 0)
        a_boot[valid_r0] = 0.5 / r0_bootstrap[valid_r0]
        epsilon_boot = np.full_like(r0_bootstrap, np.nan, dtype=float)
        epsilon_boot[valid_r0] = (self.metro.epsilon1 / self.metro.beta) * (r0_bootstrap[valid_r0] ** 2)
        length_boot = np.full_like(r0_bootstrap, np.nan, dtype=float)
        length_boot[valid_r0] = self.metro.L0 * (0.5 / r0_bootstrap[valid_r0])
        volume_boot = np.full_like(r0_bootstrap, np.nan, dtype=float)
        lattice_extent = self.metro.L0 * self.metro.L1 * self.metro.L2 * self.metro.L3
        volume_boot[valid_r0] = lattice_extent * ((0.5 / r0_bootstrap[valid_r0]) ** 4)

        np.save(self.derived_bootstrap_dir / "a.npy", a_boot)
        np.save(self.derived_bootstrap_dir / "epsilon_bar.npy", epsilon_boot)
        np.save(self.derived_bootstrap_dir / "length.npy", length_boot)
        np.save(self.derived_bootstrap_dir / "volume_r0.npy", volume_boot)

        a_mean, a_std = summarize_bootstrap(a_boot)
        eps_mean, eps_std = summarize_bootstrap(epsilon_boot)
        len_mean, len_std = summarize_bootstrap(length_boot)
        vol_mean, vol_std = summarize_bootstrap(volume_boot)

        save_json(
            summary_path,
            {
                "r0": r0_val,
                "r0_err": float(r0_err) if r0_err is not None else None,
                "a": a_mean,
                "a_err": a_std,
                "epsilon_bar": eps_mean,
                "epsilon_bar_err": eps_std,
                "length": len_mean,
                "length_err": len_std,
                "volume_r0": vol_mean,
                "volume_r0_err": vol_std,
                "bootstrap_paths": {
                    "a": str(self.derived_bootstrap_dir / "a.npy"),
                    "epsilon_bar": str(self.derived_bootstrap_dir / "epsilon_bar.npy"),
                    "length": str(self.derived_bootstrap_dir / "length.npy"),
                    "volume_r0": str(self.derived_bootstrap_dir / "volume_r0.npy"),
                },
            },
        )

        self.manifest["status"]["derived_complete"] = True
        self._save_manifest()

    def _finalize_gradient_flow(self) -> None:
        self._load_data()
        block_size = int(self.manifest.get("block_size") or 1)
        summary = run_evaluation.summarize_gradient_flow_obs(
            self.combined_gradient_flow,
            t0_target=self.gradient_flow_metadata.get("t0_target"),
            block_size=block_size,
            n_bootstrap=self.n_bootstrap,
            include_bootstrap=True,
        )
        bootstrap_samples = summary.pop("bootstrap_samples", {}) if summary else {}
        summary_path = self.gradient_flow_dir / "summary.json"
        plot_path = self.plots_dir / "gradient_flow_t2E.html"

        if summary:
            for name, samples in bootstrap_samples.items():
                np.save(self.gradient_flow_bootstrap_dir / f"{name}.npy", np.asarray(samples, dtype=float))
            summary["bootstrap_paths"] = {
                name: str(self.gradient_flow_bootstrap_dir / f"{name}.npy")
                for name in bootstrap_samples
            }
            save_json(summary_path, summary)
            save_gradient_flow_plot(plot_path, summary)
            self.manifest["gradient_flow_summary"] = str(summary_path)
            self.manifest["gradient_flow_plot"] = str(plot_path)
        else:
            save_json(summary_path, {"status": "no gradient_flow_obs.dat data available"})
            self.manifest["gradient_flow_summary"] = str(summary_path)

        self.manifest["status"]["gradient_flow_complete"] = True
        self._save_manifest()

    def _finalize_creutz(self) -> None:
        self._load_data()
        assert self.calc is not None
        summary_path = self.creutz_dir / "summary.json"
        plot_path = self.plots_dir / "creutz_ratios.html"
        diagonal_plot_path = self.plots_dir / "creutz_ratios_R_eq_T.html"

        chi: dict[str, float] = {}
        chi_err: dict[str, float] = {}
        chi_boot_paths: dict[str, str] = {}
        unique_r_set = set(self.unique_rs)
        unique_t_set = set(self.unique_ts)
        chi_pairs = [
            (float(r_val) + 0.5, float(t_val) + 0.5)
            for r_val in self.unique_rs
            if r_val + 1 in unique_r_set
            for t_val in self.unique_ts
            if t_val + 1 in unique_t_set
        ]
        status = "ok" if chi_pairs else "not enough adjacent L values for standard Creutz ratios"

        for r_mid, t_mid in chi_pairs:
            try:
                var = self.calc.get_variable("chi", R=r_mid, T=t_mid, flow_time=self.wilson_flow_time)
            except Exception:
                continue
            value = var.get()
            if value is None or not np.isfinite(value):
                continue
            key = f"{r_mid:g},{t_mid:g}"
            chi[key] = float(value)
            if var.err() is not None and np.isfinite(var.err()):
                chi_err[key] = float(var.err())
            boot = var.bootstrap()
            if boot is not None:
                filename = f"chi_R_{r_mid:g}_T_{t_mid:g}.npy".replace(".", "p")
                path = self.creutz_bootstrap_dir / filename
                np.save(path, np.asarray(boot, dtype=float))
                chi_boot_paths[key] = str(path)

        creutz_p: dict[str, float] = {}
        creutz_p_err: dict[str, float] = {}
        for r_val in [r for r in self.unique_rs if r > 0 and r % 2 == 0]:
            try:
                var = self.calc.get_variable("creutz_P", R=int(r_val), flow_time=self.wilson_flow_time)
            except Exception:
                continue
            value = var.get()
            if value is None or not np.isfinite(value):
                continue
            key = str(int(r_val))
            creutz_p[key] = float(value)
            if var.err() is not None and np.isfinite(var.err()):
                creutz_p_err[key] = float(var.err())

        payload = {
            "status": status,
            "flow_time": None if self.wilson_flow_time is None else float(self.wilson_flow_time),
            "chi": chi,
            "chi_err": chi_err,
            "chi_bootstrap_paths": chi_boot_paths,
            "creutz_P": creutz_p,
            "creutz_P_err": creutz_p_err,
        }
        save_json(summary_path, payload)
        save_creutz_plot(plot_path, payload)
        save_creutz_diagonal_plot(diagonal_plot_path, payload)
        self.manifest["creutz_summary"] = str(summary_path)
        self.manifest["creutz_plot"] = str(plot_path)
        self.manifest["creutz_diagonal_plot"] = str(diagonal_plot_path)
        self.manifest["status"]["creutz_complete"] = True
        self._save_manifest()

    def _active_v_result_rs(self) -> list[int]:
        choices_payload = self._load_v_choices()
        summary_payload = self._load_v_summary()
        active_rs: list[int] = []
        for key in choices_payload.get("choices", {}):
            r_value = int(key)
            if self._is_v_r_ignored(r_value, choices_payload):
                continue
            if summary_payload.get("results", {}).get(key) is not None:
                active_rs.append(r_value)
        return sorted(active_rs)

    def _format_r_list(self, values: list[int]) -> str:
        if not values:
            return "n/a"
        if values == list(range(values[0], values[-1] + 1)):
            return f"{values[0]}-{values[-1]}" if values[0] != values[-1] else str(values[0])
        return ", ".join(str(value) for value in values)

    def _print_result_summary(self) -> None:
        r0_result = load_json(self.r0_dir / "r0_result.json", default={}) or {}
        r0_choice = load_json(self.r0_dir / "choice.json", default={}) or {}
        derived_summary = load_json(self.derived_dir / "summary.json", default={}) or {}
        gradient_summary = load_json(self.gradient_flow_dir / "summary.json", default={}) or {}
        active_v_rs = self._active_v_result_rs()
        fit_rs = [int(value) for value in r0_result.get("fit_rs", [])]

        self.print("\nResult summary:")
        self.print(
            "  Ensemble: "
            f"beta={self.metro.beta:g}, "
            f"L={self.metro.L0}x{self.metro.L1}x{self.metro.L2}x{self.metro.L3}, "
            f"epsilon1={self.metro.epsilon1:g}, "
            f"runs={len(self.run_dirs)}"
        )
        if self.thermalization_steps is not None:
            thermalization_text = str(self.thermalization_steps)
        else:
            cuts = list(self.thermalization_steps_by_run.values())
            thermalization_text = f"per-run ({min(cuts)}-{max(cuts)})" if cuts else "per-run"
        self.print(
            f"  Analysis choices: thermalization={thermalization_text}, "
            f"wilson_flow_time={'unflowed' if self.wilson_flow_time is None else format_token(self.wilson_flow_time)}, "
            f"block_size={self.manifest.get('block_size', 'n/a')}, "
            f"bootstrap={self.n_bootstrap}"
        )
        self.print(
            f"  V(R): {len(active_v_rs)} point(s), R={self._format_r_list(active_v_rs)}"
        )
        chi2_dof = r0_result.get("chi2_dof")
        chi2_text = ""
        if chi2_dof is not None and np.isfinite(float(chi2_dof)):
            chi2_text = f", chi2/dof={float(chi2_dof):.3g}"
        self.print(
            "  r0 fit: "
            f"r_min={r0_choice.get('r_min', 'n/a')}, "
            f"fit R={self._format_r_list(fit_rs)}, "
            f"r0={format_result_value(r0_result.get('r0'), r0_result.get('r0_err'))}"
            f"{chi2_text}"
        )
        self.print(
            "  Derived: "
            f"a={format_result_value(derived_summary.get('a'), derived_summary.get('a_err'))}, "
            f"epsilon_bar={format_result_value(derived_summary.get('epsilon_bar'), derived_summary.get('epsilon_bar_err'))}, "
            f"length={format_result_value(derived_summary.get('length'), derived_summary.get('length_err'))}, "
            f"volume_r0={format_result_value(derived_summary.get('volume_r0'), derived_summary.get('volume_r0_err'))}"
        )
        t2e_0p1 = gradient_summary.get("t_over_a2_at_t2E_clover_0p1")
        t2e_0p1_fit = gradient_summary.get("t_over_a2_at_t2E_clover_0p1_weighted_fit")
        if t2e_0p1 is not None or t2e_0p1_fit is not None:
            self.print(
                "  Gradient flow: "
                "t/a^2 at t^2E_clover=0.1="
                f"{format_result_value(t2e_0p1, gradient_summary.get('t_over_a2_at_t2E_clover_0p1_err'))}, "
                "weighted four-point fit="
                f"{format_result_value(t2e_0p1_fit, gradient_summary.get('t_over_a2_at_t2E_clover_0p1_weighted_fit_err'))}"
            )

    def _write_result_summary_txt(self) -> Path:
        r0_result = load_json(self.r0_dir / "r0_result.json", default={}) or {}
        derived_summary = load_json(self.derived_dir / "summary.json", default={}) or {}
        gradient_summary = load_json(self.gradient_flow_dir / "summary.json", default={}) or {}

        t2e_0p1 = gradient_summary.get("t_over_a2_at_t2E_clover_0p1")
        t2e_0p1_err = gradient_summary.get("t_over_a2_at_t2E_clover_0p1_err")

        lines = [
            f"wilson_flow_time={'unflowed' if self.wilson_flow_time is None else format_compact_float(self.wilson_flow_time)}",
            f"r0={format_result_value_rounded_to_error(r0_result.get('r0'), r0_result.get('r0_err'))}",
            f"length={format_result_value_rounded_to_error(derived_summary.get('length'), derived_summary.get('length_err'))} fm",
            f"eps_bar={format_compact_float(derived_summary.get('epsilon_bar'))}",
            f"gf_t/a²={format_result_value_rounded_to_error(t2e_0p1, t2e_0p1_err)}",
        ]
        path = self.analysis_dir / "result_summary.txt"
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        self.manifest["result_summary_txt"] = str(path)
        self._save_manifest()
        return path

    def run(self) -> Path:
        self._finalize_wilsonloop()
        self.print("\nWilson loop analysis complete. Proceeding to finalize gradient-flow and Creutz summaries...")
        self._finalize_gradient_flow()
        self._finalize_creutz()
        self.print("\nGradient-flow and Creutz summaries complete. Proceeding to finalize r0 fit...")
        self._finalize_r0()
        self.print("\nr0 fit complete. Proceeding to compute derived quantities...")
        self._finalize_derived()
        summary_txt_path = self._write_result_summary_txt()
        self.print(f"\nFinalized analysis saved to: {self.analysis_dir}")
        self.print(f"Saved compact result summary: {summary_txt_path}")
        self._print_result_summary()
        return self.analysis_dir


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Finalize the interactive lattice analysis workflow.")
    parser.add_argument(
        "run_dirs",
        nargs="*",
        help="Exact run directories to combine, or roots containing run directories to discover.",
    )
    parser.add_argument(
        "--run-root",
        nargs="+",
        help="Discover all run directories beneath one or more roots and process them in grouped batches.",
    )
    parser.add_argument(
        "--filter",
        nargs="+",
        metavar="NAME=VALUE",
        help="Search-data-style YAML parameter filters, e.g. beta=2.4 eps1=0.0 L0=24.",
    )
    parser.add_argument(
        "--exclude",
        nargs="+",
        metavar="TEXT_OR_GLOB",
        help="Exclude run directories whose path contains this text or matches this glob.",
    )
    parser.add_argument(
        "--group-by",
        nargs="+",
        help=(
            "Optional field names used to split discovered runs, e.g. eps1 beta L0. "
            "By default, runs are grouped by compatible inputs while allowing seed, nSweep, "
            "gradient_flow.dt, and gradient_flow.t_values to vary."
        ),
    )
    parser.add_argument(
        "--require-fixed-dt",
        action="store_true",
        help="Reject or split combined runs whose GradientFlowParams.dt values differ.",
    )
    parser.add_argument(
        "--min-group-size",
        type=int,
        default=1,
        help="When grouping discovered runs, skip groups with fewer matching runs than this.",
    )
    parser.add_argument(
        "--list-groups",
        "--list-group",
        action="store_true",
        help="Print the discovered analysis groups and exit without running the final analysis.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("../data/finalized_analysis"),
        help="Root directory for finalized analysis outputs.",
    )
    parser.add_argument(
        "--plot-mode",
        default="html",
        choices=["html"],
        help="Interactive plot mode. Only HTML is currently supported.",
    )
    parser.add_argument("--load-workers", type=int, default=1, help="Parallel workers for loading input files.")
    parser.add_argument("--calc-workers", type=int, default=None, help="Reserved for future parallel calculations.")
    parser.add_argument(
        "--n-bootstrap",
        type=int,
        default=run_evaluation.DEFAULT_N_BOOTSTRAP,
        help="Number of bootstrap replicas for finalized results.",
    )
    parser.add_argument(
        "--target-force",
        type=float,
        default=run_evaluation.DEFAULT_SOMMER_TARGET,
        help="Sommer target force used in the r0 fit.",
    )
    parser.add_argument(
        "--max-r",
        type=int,
        default=None,
        help="Use only finalized V(R) points with R <= this value and ignore larger R values.",
    )
    parser.add_argument(
        "--flow-time",
        "--t-f",
        dest="wilson_flow_time",
        type=float,
        default=None,
        help=(
            "Use Wilson loops at this gradient-flow t_over_a2 value from gradient_flow_wtemp.dat "
            "for finalized V(R), r0, and Creutz analysis. Omit to use unflowed W_temp.out."
        ),
    )
    parser.add_argument(
        "--block-sizes",
        help="Comma-separated bootstrap block sizes to scan. Overrides --min-block-size/--max-block-size.",
    )
    parser.add_argument(
        "--min",
        "--min-block-size",
        dest="min_block_size",
        type=int,
        default=DEFAULT_BLOCK_SIZE_SCAN_MIN,
        help="Minimum bootstrap block size to scan when using a range.",
    )
    parser.add_argument(
        "--max",
        "--max-block-size",
        dest="max_block_size",
        type=int,
        default=DEFAULT_BLOCK_SIZE_SCAN_MAX,
        help=(
            "Maximum bootstrap block size to scan when using a range. "
            "Defaults to ceil(total post-cut configurations / 16)."
        ),
    )
    parser.add_argument(
        "--block-step",
        type=int,
        default=DEFAULT_BLOCK_SIZE_SCAN_STEP,
        help="Step for generated bootstrap block sizes.",
    )
    parser.add_argument(
        "--no-open-plots",
        action="store_false",
        dest="open_plots",
        help="Do not auto-open the HTML thermalization previews while selecting per-run cuts.",
    )
    parser.set_defaults(open_plots=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.run_root and args.run_dirs:
        parser.error("Use either positional run_dirs or --run-root roots, not both.")
    if not args.run_root and not args.run_dirs:
        parser.error("Provide run directories or use --run-root.")
    if args.min_group_size < 1:
        parser.error("--min-group-size must be at least 1.")

    try:
        block_size_scan_values = parse_bootstrap_block_sizes(args)
    except ValueError as exc:
        parser.error(str(exc))
    try:
        filter_criteria = parse_filter_tokens(args.filter)
    except ValueError as exc:
        parser.error(str(exc))

    if args.run_root:
        run_dirs = discover_run_directories_from_roots(args.run_root)
        run_dirs = filter_run_directories(run_dirs, filter_criteria)
        run_dirs = exclude_run_directories(run_dirs, args.exclude)
        if not run_dirs:
            criteria_text = " ".join(args.filter or [])
            suffix_parts = []
            if criteria_text:
                suffix_parts.append(f"matching filters: {criteria_text}")
            if args.exclude:
                suffix_parts.append(f"after exclusions: {' '.join(args.exclude)}")
            suffix = f" {'; '.join(suffix_parts)}" if suffix_parts else ""
            roots_text = ", ".join(args.run_root)
            parser.error(f"No run directories found under {roots_text}{suffix}.")
        grouped_runs = group_run_directories(
            run_dirs,
            group_by=args.group_by,
            require_fixed_dt=args.require_fixed_dt,
        )
        grouped_runs = [
            (signature, group_run_dirs)
            for signature, group_run_dirs in grouped_runs
            if len(group_run_dirs) >= args.min_group_size
        ]
        if not grouped_runs:
            parser.error(
                f"No discovered groups have at least {args.min_group_size} matching run(s)."
            )
    else:
        try:
            run_dirs, expanded_roots = resolve_selected_run_directories(args.run_dirs)
        except ValueError as exc:
            parser.error(str(exc))
        run_dirs = filter_run_directories(run_dirs, filter_criteria)
        run_dirs = exclude_run_directories(run_dirs, args.exclude)
        if not run_dirs:
            criteria_text = " ".join(args.filter or [])
            suffix_parts = []
            if criteria_text:
                suffix_parts.append(f"matching filters: {criteria_text}")
            if args.exclude:
                suffix_parts.append(f"after exclusions: {' '.join(args.exclude)}")
            suffix = f" {'; '.join(suffix_parts)}" if suffix_parts else ""
            parser.error(f"No selected run directories remain{suffix}.")
        if expanded_roots:
            grouped_runs = group_run_directories(
                run_dirs,
                group_by=args.group_by,
                require_fixed_dt=args.require_fixed_dt,
            )
            grouped_runs = [
                (signature, group_run_dirs)
                for signature, group_run_dirs in grouped_runs
                if len(group_run_dirs) >= args.min_group_size
            ]
            if not grouped_runs:
                parser.error(
                    f"No discovered groups have at least {args.min_group_size} matching run(s)."
                )
        else:
            grouped_runs = [(tuple(), run_dirs)]

    if args.list_groups:
        for index, (signature, run_dirs) in enumerate(grouped_runs, start=1):
            header = describe_group_signature(signature) if signature else "manual selection"
            print(f"[group {index}] {header}")
            for run_dir in run_dirs:
                print(f"  {run_dir}")
        return 0

    for index, (signature, run_dirs) in enumerate(grouped_runs, start=1):
        if len(grouped_runs) > 1:
            header = describe_group_signature(signature)
            print(f"\n=== Finalized analysis group {index}/{len(grouped_runs)}: {header} ===")
        runner = FinalizedAnalysisRunner(
            run_dirs=run_dirs,
            output_root=args.output_root,
            plot_mode=args.plot_mode,
            load_workers=args.load_workers,
            calc_workers=args.calc_workers,
            n_bootstrap=args.n_bootstrap,
            target_force=args.target_force,
            max_r=args.max_r,
            wilson_flow_time=args.wilson_flow_time,
            require_fixed_dt=args.require_fixed_dt,
            block_size_scan_values=block_size_scan_values,
            block_size_scan_min=args.min_block_size,
            block_size_scan_max=args.max_block_size,
            block_size_scan_step=args.block_step,
            open_plots=args.open_plots,
        )
        runner.run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
