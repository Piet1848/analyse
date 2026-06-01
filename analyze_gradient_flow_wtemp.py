#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from dataclasses import asdict
from datetime import date
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, get_origin, get_type_hints

import numpy as np

import run_evaluation
from finalized_analysis_helpers import _get_plotly, _write_figure_html
from load_input_yaml import load_gradient_flow_params, load_params


DEFAULT_N_BOOTSTRAP = 200
DEFAULT_SEED = 42
WTEMP_FILENAME = "gradient_flow_wtemp.dat"
THERMALIZATION_PREVIEW_MAX_SERIES = 12
THERMALIZATION_PREVIEW_MAX_CONFIGURATIONS = 900
BOOTSTRAP_SCAN_MAX_SERIES = 12
EXAMPLE_WILSON_PAIRS = ((1, 1), (1, 2), (2, 2), (4, 4), (8, 8), (14, 14))
EXAMPLE_WILSON_PAIR_SET = set(EXAMPLE_WILSON_PAIRS)
IGNORED_METROPOLIS_GROUP_FIELDS = {"seed", "nSweep"}
IGNORED_OUTPUT_GROUP_FIELDS = {
    "plaquette_filename",
    "W_temp_filename",
    "W_temp_L_T_pairs",
    "W_mu_nu_filename",
    "RetraceU_filename",
    "obs_filename",
    "W_mu_nu_filename",
    "t0_filename",
}
FILTER_FIELD_ALIASES = {
    "eps1": "epsilon1",
    "dt": "gradient_flow.dt",
    "gf.dt": "gradient_flow.dt",
    "gradient.dt": "gradient_flow.dt",
    "gradient_flow_dt": "gradient_flow.dt",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate mean, bootstrap error, and bootstrap samples for every "
            "Wilson loop in gradient_flow_wtemp.dat."
        )
    )
    parser.add_argument(
        "path",
        type=Path,
        help=(
            f"Path to {WTEMP_FILENAME}, a run directory containing it, "
            "or a top directory below which runs should be discovered."
        ),
    )
    parser.add_argument(
        "THERM",
        nargs="?",
        type=int,
        help=(
            "Optional thermalization cutoff. If omitted, cuts are prompted after "
            "run discovery and filtering."
        ),
    )
    parser.add_argument(
        "BLOCK",
        nargs="?",
        type=int,
        help=(
            "Optional bootstrap block size in saved measurements. If omitted, "
            "it is prompted after thermalization selection."
        ),
    )
    parser.add_argument(
        "--n-bootstrap",
        "--bootstrap-samples",
        dest="n_bootstrap",
        type=int,
        default=DEFAULT_N_BOOTSTRAP,
        help=f"Number of bootstrap replicas to generate. Default: {DEFAULT_N_BOOTSTRAP}.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_SEED,
        help=f"Random seed for bootstrap resampling. Default: {DEFAULT_SEED}.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        help=(
            "Output directory. Default: <input directory>/gradient_flow_wtemp_analysis. "
            "For discovered groups, one subdirectory is written per parameter group."
        ),
    )
    parser.add_argument(
        "--date-min",
        help=(
            "Minimum run date to include, read from path components. "
            "Formats: YYYY-MM-DD, YYYY_MM_DD, or YYYYMMDD; suffixes like 20260523_a100 are accepted."
        ),
    )
    parser.add_argument(
        "--date-max",
        help=(
            "Maximum run date to include, read from path components. "
            "Formats: YYYY-MM-DD, YYYY_MM_DD, or YYYYMMDD; suffixes like 20260523_a100 are accepted."
        ),
    )
    parser.add_argument(
        "--date-range",
        nargs=2,
        metavar=("DATE_MIN", "DATE_MAX"),
        help=(
            "Inclusive run-date range to include. Formats: YYYY-MM-DD, YYYY_MM_DD, or YYYYMMDD."
        ),
    )
    parser.add_argument(
        "--min-group-size",
        type=int,
        default=1,
        help="Skip discovered parameter groups with fewer runs than this. Default: 1.",
    )
    parser.add_argument(
        "--list-groups",
        action="store_true",
        help="Print discovered parameter groups and exit without writing analysis output.",
    )
    parser.add_argument(
        "--filter",
        nargs="+",
        metavar="NAME=VALUE",
        help="YAML parameter filters applied before grouping, e.g. beta=2.4 eps1=0.0 L0=24.",
    )
    return parser.parse_args()


def build_filter_field_map() -> dict[str, tuple[str, Any]]:
    from load_input_yaml import GaugeObservableParams, GradientFlowParams, MetropolisParams

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


def run_matches_filter(run_dir: Path, criteria: dict[str, Any]) -> bool:
    if not criteria:
        return True
    metro, gauge = load_params(str(run_dir / "input.yaml"))
    gradient_flow = load_gradient_flow_params(str(run_dir / "input.yaml"))
    for name, expected in criteria.items():
        block_name, _ = FILTER_FIELD_MAP[name]
        if block_name == "metro":
            obj = metro
            attr_name = name
        elif block_name == "gauge":
            obj = gauge
            attr_name = name
        else:
            obj = gradient_flow
            attr_name = name.split(".", 1)[1] if "." in name else name
        actual = getattr(obj, attr_name)
        if isinstance(actual, (list, dict)) and isinstance(expected, str):
            if expected not in str(actual):
                return False
            continue
        if actual != expected:
            return False
    return True


def normalize_path(path: Path) -> Path:
    return path.expanduser().resolve()


def parse_path_date(value: str) -> date:
    cleaned = value.strip()
    for fmt in ("%Y-%m-%d", "%Y_%m_%d", "%Y%m%d"):
        try:
            return datetime.strptime(cleaned, fmt).date()
        except ValueError:
            pass
    raise ValueError(f"invalid date {value!r}; expected YYYY-MM-DD or YYYYMMDD")


def path_dates(path: Path) -> list[date]:
    dates: list[date] = []
    seen: set[date] = set()
    pattern = re.compile(r"(?<!\d)(20\d{2})[-_]?([01]\d)[-_]?([0-3]\d)(?!\d)")
    for part in path.parts:
        for match in pattern.finditer(part):
            try:
                current = date(int(match.group(1)), int(match.group(2)), int(match.group(3)))
            except ValueError:
                continue
            if current not in seen:
                dates.append(current)
                seen.add(current)
    return dates


def path_matches_date_range(
    path: Path,
    date_min: date | None,
    date_max: date | None,
) -> bool:
    if date_min is None and date_max is None:
        return True
    for current in path_dates(path):
        if date_min is not None and current < date_min:
            continue
        if date_max is not None and current > date_max:
            continue
        return True
    return False


def discover_run_dirs(top_path: Path) -> list[Path]:
    top_path = normalize_path(top_path)
    if not top_path.is_dir():
        raise ValueError(f"top path is not a directory: {top_path}")
    return sorted({path.parent for path in top_path.rglob("input.yaml") if path.is_file()})


def gradient_wtemp_path_for_run(run_dir: Path) -> Path:
    yaml_path = run_dir / "input.yaml"
    try:
        gradient_flow = load_gradient_flow_params(str(yaml_path))
        filename = gradient_flow.W_temp_filename or WTEMP_FILENAME
    except Exception:
        filename = WTEMP_FILENAME
    return run_dir / filename


def jsonable_value(value: Any) -> Any:
    if isinstance(value, tuple):
        return [jsonable_value(item) for item in value]
    if isinstance(value, list):
        return [jsonable_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): jsonable_value(val) for key, val in sorted(value.items())}
    return value


def freeze_group_value(value: Any) -> Any:
    if isinstance(value, list):
        return tuple(freeze_group_value(item) for item in value)
    if isinstance(value, tuple):
        return tuple(freeze_group_value(item) for item in value)
    if isinstance(value, dict):
        return tuple(
            (str(key), freeze_group_value(val))
            for key, val in sorted(value.items())
        )
    return value


def group_signature_for_run(run_dir: Path) -> tuple[tuple[str, Any], ...]:
    metro, gauge = load_params(str(run_dir / "input.yaml"))
    gradient_flow = load_gradient_flow_params(str(run_dir / "input.yaml"))

    fields: dict[str, Any] = {}
    for key, value in asdict(metro).items():
        if key not in IGNORED_METROPOLIS_GROUP_FIELDS:
            fields[f"MetropolisParams.{key}"] = freeze_group_value(value)
    for key, value in asdict(gauge).items():
        if key not in IGNORED_OUTPUT_GROUP_FIELDS:
            fields[f"GaugeObservableParams.{key}"] = freeze_group_value(value)
    for key, value in asdict(gradient_flow).items():
        if key not in IGNORED_OUTPUT_GROUP_FIELDS:
            fields[f"GradientFlowParams.{key}"] = freeze_group_value(value)

    return tuple(sorted(fields.items()))


def short_group_token(signature: tuple[tuple[str, Any], ...]) -> str:
    selected: list[str] = []
    lookup = dict(signature)
    for key in (
        "MetropolisParams.beta",
        "MetropolisParams.L0",
        "MetropolisParams.epsilon1",
        "GradientFlowParams.dt",
    ):
        if key in lookup:
            selected.append(f"{key.rsplit('.', 1)[1]}_{filename_token(lookup[key])}")
    text = "__".join(selected) if selected else "manual"
    digest_source = json.dumps(signature, sort_keys=True, default=str).encode("utf-8")
    digest = hashlib.sha1(digest_source).hexdigest()[:8]
    safe_text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_") or "group"
    return f"{safe_text}__{digest}"


def describe_group(signature: tuple[tuple[str, Any], ...]) -> str:
    lookup = dict(signature)
    parts = []
    for key in (
        "MetropolisParams.beta",
        "MetropolisParams.L0",
        "MetropolisParams.L1",
        "MetropolisParams.L2",
        "MetropolisParams.L3",
        "MetropolisParams.epsilon1",
        "GradientFlowParams.dt",
    ):
        if key in lookup:
            parts.append(f"{key.rsplit('.', 1)[1]}={lookup[key]}")
    return ", ".join(parts) if parts else "manual"


def differing_group_fields(groups: list[dict[str, Any]]) -> list[str]:
    if len(groups) <= 1:
        return []
    signatures = [dict(group["signature"]) for group in groups]
    field_names = sorted({name for signature in signatures for name in signature})
    differing: list[str] = []
    for name in field_names:
        values = {repr(signature.get(name)) for signature in signatures}
        if len(values) > 1:
            differing.append(name)
    return differing


def format_group_difference(signature: tuple[tuple[str, Any], ...], fields: list[str]) -> str:
    if not fields:
        return ""
    lookup = dict(signature)
    parts = [f"{field}={format_difference_value(lookup.get(field))}" for field in fields]
    return "; ".join(parts)


def compact_int_range(values: list[int]) -> str:
    unique = sorted(set(int(value) for value in values))
    if not unique:
        return "[]"
    if unique == list(range(unique[0], unique[-1] + 1)):
        return f"{unique[0]}:{unique[-1]}"
    return ",".join(str(value) for value in unique)


def format_difference_value(value: Any) -> str:
    jsonable = jsonable_value(value)
    if (
        isinstance(jsonable, list)
        and jsonable
        and all(isinstance(item, list) and len(item) == 2 for item in jsonable)
    ):
        try:
            first = [int(item[0]) for item in jsonable]
            second = [int(item[1]) for item in jsonable]
        except (TypeError, ValueError):
            return str(jsonable)
        return f"{compact_int_range(first)} x {compact_int_range(second)} ({len(jsonable)} pairs)"
    return str(jsonable)


def pair_coverage_for_run(run_dir: Path) -> dict[str, Any]:
    _metro, gauge = load_params(str(run_dir / "input.yaml"))
    pairs = [(int(r_val), int(t_val)) for r_val, t_val in gauge.W_temp_L_T_pairs]
    if not pairs:
        return {"pairs": [], "label": "none", "n_pairs": 0}
    rs = [pair[0] for pair in pairs]
    ts = [pair[1] for pair in pairs]
    return {
        "pairs": [[r_val, t_val] for r_val, t_val in pairs],
        "label": f"{compact_int_range(rs)} x {compact_int_range(ts)} ({len(pairs)} pairs)",
        "n_pairs": len(pairs),
    }


def pair_coverage_summary(run_dirs: list[Path]) -> dict[str, Any]:
    by_run: dict[str, dict[str, Any]] = {}
    counts: dict[str, int] = {}
    for run_dir in run_dirs:
        coverage = pair_coverage_for_run(run_dir)
        by_run[str(run_dir)] = coverage
        counts[coverage["label"]] = counts.get(coverage["label"], 0) + 1
    return {
        "by_run": by_run,
        "groups": [
            {"label": label, "n_runs": count}
            for label, count in sorted(counts.items())
        ],
    }


def print_pair_coverage_if_mixed(group: dict[str, Any]) -> None:
    run_dirs = group.get("run_dirs") or []
    if not run_dirs:
        return
    summary = pair_coverage_summary(run_dirs)
    if len(summary["groups"]) <= 1:
        return
    print("  pair coverage:")
    for row in summary["groups"]:
        print(f"    {row['label']}: {row['n_runs']} run(s)")


def resolve_date_filters(args: argparse.Namespace) -> tuple[date | None, date | None]:
    date_min = args.date_min
    date_max = args.date_max
    if args.date_range:
        if date_min is not None or date_max is not None:
            raise ValueError("Use either --date-range or --date-min/--date-max, not both.")
        date_min, date_max = args.date_range

    parsed_min = parse_path_date(date_min) if date_min else None
    parsed_max = parse_path_date(date_max) if date_max else None
    if parsed_min is not None and parsed_max is not None and parsed_min > parsed_max:
        raise ValueError("--date-min must be <= --date-max")
    return parsed_min, parsed_max


def resolve_input_groups(
    path: Path,
    date_min: date | None,
    date_max: date | None,
    min_group_size: int,
    filter_criteria: dict[str, Any],
) -> tuple[bool, list[dict[str, Any]]]:
    path = normalize_path(path)
    if path.is_file():
        if filter_criteria:
            raise ValueError("--filter requires a run directory or top path with input.yaml files")
        return False, [
            {
                "signature": tuple(),
                "run_dirs": [],
                "input_paths": [path],
            }
        ]

    if not path.is_dir():
        raise ValueError(f"input path does not exist: {path}")

    if (path / "input.yaml").is_file():
        run_dirs = [path]
        discovered = False
    else:
        run_dirs = discover_run_dirs(path)
        discovered = True

    run_dirs = [
        run_dir
        for run_dir in run_dirs
        if path_matches_date_range(run_dir, date_min, date_max)
    ]
    filtered_run_dirs: list[Path] = []
    for run_dir in run_dirs:
        try:
            if run_matches_filter(run_dir, filter_criteria):
                filtered_run_dirs.append(run_dir)
        except Exception as exc:
            raise ValueError(f"failed to apply filters to {run_dir / 'input.yaml'}: {exc}") from exc
    run_dirs = filtered_run_dirs
    if not run_dirs:
        raise ValueError("no run directories remain after date and parameter filtering")

    grouped: dict[tuple[tuple[str, Any], ...], list[Path]] = {}
    for run_dir in run_dirs:
        try:
            signature = group_signature_for_run(run_dir)
        except Exception as exc:
            raise ValueError(f"failed to read parameters from {run_dir / 'input.yaml'}: {exc}") from exc
        grouped.setdefault(signature, []).append(run_dir)

    groups: list[dict[str, Any]] = []
    for signature, group_run_dirs in sorted(grouped.items(), key=lambda item: describe_group(item[0])):
        if len(group_run_dirs) < min_group_size:
            continue
        input_paths = [gradient_wtemp_path_for_run(run_dir) for run_dir in group_run_dirs]
        missing = [str(input_path) for input_path in input_paths if not input_path.exists()]
        if missing:
            raise ValueError(
                "missing gradient-flow Wilson-loop file(s): " + ", ".join(missing)
            )
        groups.append(
            {
                "signature": signature,
                "run_dirs": group_run_dirs,
                "input_paths": input_paths,
            }
        )
    if not groups:
        raise ValueError(f"no parameter groups have at least {min_group_size} run(s)")
    return discovered, groups


def split_table_line(line: str, delimiter: str | None) -> list[str]:
    if delimiter == ",":
        return [token.strip() for token in next(csv.reader([line], skipinitialspace=True))]
    return line.split()


def looks_like_column_name(token: str) -> bool:
    return any(ch.isalpha() for ch in token) and all(
        ch.isalnum() or ch in {"_", "/", "#"} for ch in token
    )


def parse_integer_token(token: str, path: Path, line_number: int, name: str) -> int:
    try:
        value = float(token)
    except ValueError as exc:
        raise ValueError(f"{path}:{line_number}: could not parse {name}") from exc

    integer_value = int(value)
    if value != integer_value:
        raise ValueError(f"{path}:{line_number}: {name} is not an integer")
    return integer_value


def canonical_header(tokens: list[str] | None, width: int) -> list[str]:
    if tokens is not None and len(tokens) == width:
        return tokens
    if width == 5:
        return ["conf_id", "t_over_a2", "L", "T", "W_temp"]
    raise ValueError(
        "expected gradient-flow W_temp rows with 5 columns: "
        "conf_id t_over_a2 L T W_temp"
    )


def iter_gradient_flow_wtemp_rows(path: Path):
    delimiter: str | None = None
    header_tokens: list[str] | None = None
    header: list[str] | None = None

    with path.open("r", encoding="utf-8-sig") as stream:
        for line_number, line in enumerate(stream, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                tokens = stripped[1:].strip().split()
                if tokens and all(looks_like_column_name(token) for token in tokens):
                    header_tokens = tokens
                continue

            if delimiter is None:
                delimiter = "," if "," in stripped else None
            fields = split_table_line(stripped, delimiter)
            if not fields:
                continue

            if header is None and all(looks_like_column_name(token) for token in fields):
                header_tokens = fields
                continue

            if header is None:
                header = canonical_header(header_tokens, len(fields))
            if len(fields) != len(header):
                raise ValueError(
                    f"{path}:{line_number}: expected {len(header)} columns, got {len(fields)}"
                )

            row = dict(zip(header, fields))
            try:
                conf_id = float(row["conf_id"])
                flow_time = float(row["t_over_a2"])
                wilson_loop = float(row["W_temp"])
            except (KeyError, ValueError) as exc:
                raise ValueError(
                    f"{path}:{line_number}: could not parse gradient-flow W_temp row"
                ) from exc

            yield (
                conf_id,
                flow_time,
                parse_integer_token(row["L"], path, line_number, "L"),
                parse_integer_token(row["T"], path, line_number, "T"),
                wilson_loop,
            )


def infer_rows_per_configuration(path: Path) -> int:
    first_conf_id: float | None = None
    count = 0
    for conf_id, _flow_time, _l_size, _t_size, _wilson_loop in iter_gradient_flow_wtemp_rows(path):
        if first_conf_id is None:
            first_conf_id = conf_id
        elif conf_id != first_conf_id:
            return count
        count += 1
    return count


def read_gradient_flow_wtemp(
    path: Path,
    therm: int,
    key_filter=None,
) -> tuple[dict[tuple[float, int, int], np.ndarray], dict[tuple[float, int, int], tuple[float, float]]]:
    table = np.loadtxt(path, comments="#", dtype=np.float64)
    if table.size == 0:
        return {}, {}
    if table.ndim == 1:
        table = table.reshape(1, -1)
    if table.shape[1] != 5:
        raise ValueError(
            f"{path}: expected 5 columns: conf_id t_over_a2 L T W_temp"
        )

    conf_ids = table[:, 0]
    change_idx = np.flatnonzero(conf_ids != conf_ids[0])
    rows_per_conf = int(change_idx[0]) if change_idx.size else int(len(conf_ids))
    if rows_per_conf <= 0:
        return {}, {}

    complete_rows = (len(conf_ids) // rows_per_conf) * rows_per_conf
    if complete_rows <= 0:
        return {}, {}
    if complete_rows != len(conf_ids):
        table = table[:complete_rows]
        conf_ids = conf_ids[:complete_rows]

    n_configurations = complete_rows // rows_per_conf
    first_block = table[:rows_per_conf]
    key_order = [
        (float(row[1]), int(row[2]), int(row[3]))
        for row in first_block
    ]
    if len(set(key_order)) != len(key_order):
        raise ValueError(f"{path}: duplicate (t_over_a2, L, T) rows in first configuration")

    conf_by_configuration = conf_ids.reshape(n_configurations, rows_per_conf)[:, 0]
    keep_configurations = conf_by_configuration >= float(therm)
    if not np.any(keep_configurations):
        return {}, {}

    w_matrix = table[:, 4].reshape(n_configurations, rows_per_conf)
    kept_conf_ids = conf_by_configuration[keep_configurations]
    conf_range = (float(np.min(kept_conf_ids)), float(np.max(kept_conf_ids)))

    values_by_loop: dict[tuple[float, int, int], np.ndarray] = {}
    conf_range_by_loop: dict[tuple[float, int, int], tuple[float, float]] = {}
    selected_columns = [
        (idx, key)
        for idx, key in enumerate(key_order)
        if key_filter is None or key_filter(key)
    ]
    for idx, key in selected_columns:
        values_by_loop[key] = np.asarray(w_matrix[keep_configurations, idx], dtype=float)
        conf_range_by_loop[key] = conf_range

    return values_by_loop, conf_range_by_loop


def read_combined_gradient_flow_wtemp(
    input_paths: list[Path],
    thermalization_by_input: dict[Path, int],
    *,
    label: str | None = None,
    key_filter=None,
) -> tuple[dict[tuple[float, int, int], np.ndarray], dict[tuple[float, int, int], tuple[float, float]]]:
    chunks_by_loop: dict[tuple[float, int, int], list[np.ndarray]] = {}
    conf_range_by_loop: dict[tuple[float, int, int], tuple[float, float]] = {}
    for index, input_path in enumerate(input_paths, start=1):
        if label:
            print(f"  {label}: loading file {index}/{len(input_paths)}: {input_path}")
        values_by_loop, current_conf_range_by_loop = read_gradient_flow_wtemp(
            input_path,
            thermalization_by_input[input_path],
            key_filter=key_filter,
        )
        for loop_key, values in values_by_loop.items():
            chunks_by_loop.setdefault(loop_key, []).append(np.asarray(values, dtype=float))
            first_conf, last_conf = current_conf_range_by_loop[loop_key]
            if loop_key not in conf_range_by_loop:
                conf_range_by_loop[loop_key] = (first_conf, last_conf)
            else:
                old_first, old_last = conf_range_by_loop[loop_key]
                conf_range_by_loop[loop_key] = (
                    min(old_first, first_conf),
                    max(old_last, last_conf),
                )
    combined_values = {
        loop_key: np.concatenate(chunks).astype(float, copy=False)
        for loop_key, chunks in chunks_by_loop.items()
        if chunks
    }
    return combined_values, conf_range_by_loop


def build_bootstrap_plan(
    n_samples: int,
    block_size: int,
    n_bootstrap: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int, int]:
    starts = np.arange(0, n_samples, block_size, dtype=np.int64)
    lengths = np.minimum(block_size, n_samples - starts).astype(np.float64)
    n_blocks = int(len(starts))
    sampled_blocks = rng.integers(0, n_blocks, size=(n_bootstrap, n_blocks))
    tail = int(n_samples - ((n_samples // block_size) * block_size))
    return starts, lengths, sampled_blocks, n_blocks, tail


def bootstrap_means(
    values: np.ndarray,
    plan: tuple[np.ndarray, np.ndarray, np.ndarray, int, int],
) -> tuple[np.ndarray, int, int]:
    starts, lengths, sampled_blocks, n_blocks, tail = plan
    block_sums = np.add.reduceat(values.astype(np.float64, copy=False), starts)
    samples = (
        block_sums[sampled_blocks].sum(axis=1)
        / lengths[sampled_blocks].sum(axis=1)
    )
    return samples.astype(float, copy=False), n_blocks, tail


def filename_token(value: float | int | str) -> str:
    text = f"{value:.12g}" if isinstance(value, float) else str(value)
    return text.replace("-", "m").replace(".", "p").replace("+", "")


def relative_to_output(path: Path, output_dir: Path) -> str:
    try:
        return str(path.relative_to(output_dir))
    except ValueError:
        return str(path)


def prompt_thermalization_for_group(
    group: dict[str, Any],
    default_cut: int,
) -> dict[Path, int]:
    run_dirs: list[Path] = group["run_dirs"]
    input_paths: list[Path] = group["input_paths"]
    if not run_dirs:
        input_path = input_paths[0]
        while True:
            response = input(
                f"Thermalization cut for {input_path.name} is {default_cut}. "
                "Press Enter to accept, or type a new integer: "
            ).strip()
            if response == "" or response.lower() in {"y", "yes"}:
                return {input_path: int(default_cut)}
            try:
                value = int(response)
            except ValueError:
                print(f"Invalid thermalization cut: {response!r}")
                continue
            if value < 0:
                print("Thermalization cut must be non-negative.")
                continue
            return {input_path: value}

    cuts_by_input: dict[Path, int] = {}
    print("Selected run directories:")
    for index, run_dir in enumerate(run_dirs, start=1):
        print(f"  {index}. {run_dir}")

    for index, (run_dir, input_path) in enumerate(zip(run_dirs, input_paths, strict=True), start=1):
        while True:
            response = input(
                f"Thermalization cut for {run_dir.name} is {default_cut}. "
                "Press Enter to accept, type a new integer, or append * to apply to all runs: "
            ).strip()
            if response == "" or response.lower() in {"y", "yes"}:
                cuts_by_input[input_path] = int(default_cut)
                break
            apply_to_all = response.endswith("*")
            if apply_to_all:
                response = response[:-1].strip()
            try:
                value = int(response)
            except ValueError:
                print(f"Invalid thermalization cut: {response!r}")
                continue
            if value < 0:
                print("Thermalization cut must be non-negative.")
                continue
            if apply_to_all:
                return {path: value for path in input_paths}
            cuts_by_input[input_path] = value
            break
        print(f"  Run {index}/{len(run_dirs)}: {run_dir.name} -> {cuts_by_input[input_path]}")
    return cuts_by_input


def prompt_bootstrap_block_size(default_block_size: int) -> int:
    while True:
        response = input(
            "Bootstrap block size in saved configurations "
            f"is {int(default_block_size)}. Press Enter to accept, or type a new integer: "
        ).strip()
        if response == "" or response.lower() in {"y", "yes"}:
            return int(default_block_size)
        try:
            value = int(response)
        except ValueError:
            print(f"Invalid bootstrap block size: {response!r}")
            continue
        if value < 1:
            print("Bootstrap block size must be at least 1.")
            continue
        return value


def is_example_loop_pair(l_size: int, t_size: int) -> bool:
    return (int(l_size), int(t_size)) in EXAMPLE_WILSON_PAIR_SET


def select_example_loop_keys(
    available_keys,
    max_series: int,
) -> list[tuple[float, int, int]]:
    key_set = {
        (float(flow_time), int(l_size), int(t_size))
        for flow_time, l_size, t_size in available_keys
    }
    flow_times = sorted({key[0] for key in key_set})
    selected: list[tuple[float, int, int]] = []
    for flow_time in flow_times:
        for l_size, t_size in EXAMPLE_WILSON_PAIRS:
            key = (float(flow_time), int(l_size), int(t_size))
            if key in key_set:
                selected.append(key)
                if len(selected) >= max_series:
                    return selected

    for key in sorted(key_set):
        if key not in selected:
            selected.append(key)
            if len(selected) >= max_series:
                return selected
    return selected


def save_decision_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True, default=str)
        stream.write("\n")


def save_thermalization_preview_plot(
    path: Path,
    group: dict[str, Any],
    suggested_cut: int,
) -> None:
    go = _get_plotly()
    from plotly.subplots import make_subplots

    series_by_key: dict[tuple[float, int, int], list[dict[str, Any]]] = {}
    run_dirs: list[Path] = group["run_dirs"]
    input_paths: list[Path] = group["input_paths"]
    selected_keys: list[tuple[float, int, int]] | None = None

    for index, input_path in enumerate(input_paths):
        run_label = run_dirs[index].name if index < len(run_dirs) else input_path.parent.name
        rows_per_conf = infer_rows_per_configuration(input_path)
        if rows_per_conf <= 0:
            continue
        max_rows = rows_per_conf * THERMALIZATION_PREVIEW_MAX_CONFIGURATIONS
        print(
            "  Building thermalization preview from "
            f"{input_path} using up to {THERMALIZATION_PREVIEW_MAX_CONFIGURATIONS} saved configurations..."
        )
        values_by_key: dict[tuple[float, int, int], list[float]] = {}
        steps_by_key: dict[tuple[float, int, int], list[float]] = {}
        preview_table = np.loadtxt(
            input_path,
            comments="#",
            dtype=np.float64,
            max_rows=max_rows,
        )
        if preview_table.ndim == 1:
            preview_table = preview_table.reshape(1, -1)
        for conf_id, flow_time, l_size_raw, t_size_raw, wilson_loop in preview_table:
            l_size = int(l_size_raw)
            t_size = int(t_size_raw)
            key = (flow_time, l_size, t_size)
            if not is_example_loop_pair(l_size, t_size):
                continue
            if selected_keys is not None and key not in selected_keys:
                continue
            steps_by_key.setdefault(key, []).append(conf_id)
            values_by_key.setdefault(key, []).append(wilson_loop)
            if selected_keys is None and len(values_by_key) >= THERMALIZATION_PREVIEW_MAX_SERIES:
                selected_keys = select_example_loop_keys(
                    values_by_key,
                    THERMALIZATION_PREVIEW_MAX_SERIES,
                )
        for key, values in values_by_key.items():
            if selected_keys is not None and key not in selected_keys:
                continue
            arr = np.asarray(values, dtype=float)
            if arr.size == 0:
                continue
            running = np.cumsum(arr) / np.arange(1, arr.size + 1, dtype=float)
            series_by_key.setdefault(key, []).append(
                {
                    "run_label": run_label,
                    "steps": steps_by_key[key],
                    "values": values,
                    "running_mean": running.tolist(),
                }
            )

    if not series_by_key:
        raise ValueError("no thermalization preview data available")

    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.09,
        subplot_titles=("Raw W(L,T) vs configuration", "Running mean of W(L,T)"),
    )
    trace_indices_by_key: dict[tuple[float, int, int], list[int]] = {}
    ordered_keys = sorted(series_by_key)
    for key_index, key in enumerate(ordered_keys):
        visible = key_index == 0
        trace_indices_by_key[key] = []
        for row in series_by_key[key]:
            steps = [float(step) for step in row["steps"]]
            values = [float(value) for value in row["values"]]
            running = [float(value) for value in row["running_mean"]]
            name = str(row["run_label"])
            trace_indices_by_key[key].append(len(fig.data))
            fig.add_trace(
                go.Scatter(
                    x=steps,
                    y=values,
                    mode="lines+markers",
                    name=name,
                    visible=visible,
                    opacity=0.55,
                    legendgroup=name,
                ),
                row=1,
                col=1,
            )
            trace_indices_by_key[key].append(len(fig.data))
            fig.add_trace(
                go.Scatter(
                    x=steps,
                    y=running,
                    mode="lines",
                    name=f"{name} running",
                    visible=visible,
                    legendgroup=name,
                    showlegend=False,
                ),
                row=2,
                col=1,
            )

    shapes = [
        {
            "type": "line",
            "xref": "x",
            "yref": "paper",
            "x0": int(suggested_cut),
            "x1": int(suggested_cut),
            "y0": 0,
            "y1": 1,
            "line": {"color": "red", "width": 1.5, "dash": "dash"},
        }
    ]
    buttons = []
    for key in ordered_keys:
        visible = [False] * len(fig.data)
        for trace_index in trace_indices_by_key[key]:
            visible[trace_index] = True
        flow_time, l_size, t_size = key
        buttons.append(
            {
                "label": f"t={flow_time:g}, L={l_size}, T={t_size}",
                "method": "update",
                "args": [
                    {"visible": visible},
                    {"title": f"Thermalization preview: t/a^2={flow_time:g}, L={l_size}, T={t_size}"},
                ],
            }
        )

    first_key = ordered_keys[0]
    fig.update_layout(
        title=f"Thermalization preview: t/a^2={first_key[0]:g}, L={first_key[1]}, T={first_key[2]}",
        template="plotly_white",
        xaxis2_title="configuration id",
        yaxis_title="W(L,T)",
        yaxis2_title="running mean",
        shapes=shapes,
        updatemenus=[
            {
                "buttons": buttons,
                "direction": "down",
                "x": 1.0,
                "xanchor": "right",
                "y": 1.16,
                "yanchor": "top",
            }
        ],
        annotations=[
            {
                "text": f"suggested thermalization cut = {int(suggested_cut)}",
                "xref": "paper",
                "yref": "paper",
                "x": 0.01,
                "y": 1.10,
                "showarrow": False,
            }
        ],
    )
    _write_figure_html(fig, path)


def block_size_scan_values(max_samples: int, default_block_size: int) -> list[int]:
    candidates = [
        *range(1, 11),
        *range(20, 501, 10),
        *range(525, 2001, 25),
        int(default_block_size),
        int(default_block_size) * 2,
    ]
    if max_samples > 0:
        candidates.append(max_samples)
    return sorted({value for value in candidates if value >= 1 and (max_samples <= 0 or value <= max_samples)})


def build_bootstrap_block_size_scan_records(
    input_paths: list[Path],
    thermalization_by_input: dict[Path, int],
    block_sizes: list[int],
    n_bootstrap: int,
    seed: int,
) -> list[dict[str, Any]]:
    values_by_loop, _conf_range_by_loop = read_combined_gradient_flow_wtemp(
        input_paths,
        thermalization_by_input,
        label="block-size scan",
        key_filter=lambda key: is_example_loop_pair(key[1], key[2]),
    )
    selected_keys = select_example_loop_keys(values_by_loop, BOOTSTRAP_SCAN_MAX_SERIES)
    records: list[dict[str, Any]] = []
    for block_size in block_sizes:
        rng = np.random.default_rng(seed)
        plans: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray, int, int]] = {}
        for flow_time, l_size, t_size in selected_keys:
            values = np.asarray(values_by_loop[(flow_time, l_size, t_size)], dtype=float)
            if values.size == 0:
                continue
            plan = plans.get(int(values.size))
            if plan is None:
                plan = build_bootstrap_plan(int(values.size), int(block_size), n_bootstrap, rng)
                plans[int(values.size)] = plan
            samples, n_blocks, _tail = bootstrap_means(values, plan)
            records.append(
                {
                    "t_over_a2": float(flow_time),
                    "L": int(l_size),
                    "T": int(t_size),
                    "n_measurements": int(values.size),
                    "block_size": int(block_size),
                    "n_blocks": int(n_blocks),
                    "mean": float(np.mean(values, dtype=np.float64)),
                    "error": float(np.std(samples)),
                }
            )
    return records


def save_bootstrap_block_size_scan_plot(
    path: Path,
    records: list[dict[str, Any]],
    recommended_block_size: int,
) -> None:
    go = _get_plotly()
    fig = go.Figure()
    grouped: dict[tuple[float, int, int], list[dict[str, Any]]] = {}
    for record in records:
        key = (float(record["t_over_a2"]), int(record["L"]), int(record["T"]))
        grouped.setdefault(key, []).append(record)

    for flow_time, l_size, t_size in sorted(grouped):
        rows = sorted(grouped[(flow_time, l_size, t_size)], key=lambda row: int(row["block_size"]))
        fig.add_trace(
            go.Scatter(
                x=[int(row["block_size"]) for row in rows],
                y=[float(row["error"]) for row in rows],
                mode="lines+markers",
                name=f"t={flow_time:g}, L={l_size}, T={t_size}",
            )
        )

    fig.update_layout(
        title="Bootstrap error vs block size",
        template="plotly_white",
        xaxis_title="Bootstrap block size",
        yaxis_title="Bootstrap error of W(L,T)",
        shapes=[
            {
                "type": "line",
                "xref": "x",
                "yref": "paper",
                "x0": int(recommended_block_size),
                "x1": int(recommended_block_size),
                "y0": 0,
                "y1": 1,
                "line": {"color": "red", "width": 1.5, "dash": "dash"},
            }
        ],
        annotations=[
            {
                "text": f"default block size = {int(recommended_block_size)}",
                "xref": "paper",
                "yref": "paper",
                "x": 0.01,
                "y": 1.08,
                "showarrow": False,
            }
        ],
    )
    _write_figure_html(fig, path)


def save_wilson_loop_plot(path: Path, records: list[dict[str, Any]]) -> None:
    go = _get_plotly()
    fig = go.Figure()

    by_flow: dict[float, list[dict[str, Any]]] = {}
    for record in records:
        by_flow.setdefault(float(record["t_over_a2"]), []).append(record)

    trace_indices_by_flow: dict[float, list[int]] = {}
    for flow_index, flow_time in enumerate(sorted(by_flow)):
        visible = flow_index == 0
        trace_indices_by_flow[flow_time] = []
        rows_by_l: dict[int, list[dict[str, Any]]] = {}
        for record in by_flow[flow_time]:
            rows_by_l.setdefault(int(record["L"]), []).append(record)

        for l_size in sorted(rows_by_l):
            rows = sorted(rows_by_l[l_size], key=lambda item: int(item["T"]))
            trace_indices_by_flow[flow_time].append(len(fig.data))
            fig.add_trace(
                go.Scatter(
                    x=[int(row["T"]) for row in rows],
                    y=[float(row["mean"]) for row in rows],
                    error_y={
                        "type": "data",
                        "array": [float(row["error"]) for row in rows],
                        "visible": True,
                    },
                    customdata=[
                        [
                            int(row["L"]),
                            int(row["T"]),
                            int(row["n_measurements"]),
                            float(row["error"]),
                        ]
                        for row in rows
                    ],
                    hovertemplate=(
                        "L=%{customdata[0]}<br>"
                        "T=%{customdata[1]}<br>"
                        "n=%{customdata[2]}<br>"
                        "mean=%{y:g}<br>"
                        "bootstrap error=%{customdata[3]:g}"
                        "<extra>t/a^2 %{fullData.name}</extra>"
                    ),
                    mode="lines+markers",
                    name=f"L={l_size}",
                    legendgroup=f"flow_{flow_time:g}",
                    visible=visible,
                )
            )

    buttons = []
    for flow_time in sorted(by_flow):
        visible = [False] * len(fig.data)
        for trace_index in trace_indices_by_flow[flow_time]:
            visible[trace_index] = True
        buttons.append(
            {
                "label": f"t/a^2={flow_time:g}",
                "method": "update",
                "args": [
                    {"visible": visible},
                    {"title": f"Gradient-flow Wilson loops at t/a^2={flow_time:g}"},
                ],
            }
        )

    if buttons:
        fig.update_layout(
            updatemenus=[
                {
                    "buttons": buttons,
                    "direction": "down",
                    "x": 1.0,
                    "xanchor": "right",
                    "y": 1.16,
                    "yanchor": "top",
                }
            ]
        )

    first_flow = min(by_flow) if by_flow else None
    title = (
        f"Gradient-flow Wilson loops at t/a^2={first_flow:g}"
        if first_flow is not None
        else "Gradient-flow Wilson loops"
    )
    fig.update_layout(
        title=title,
        template="plotly_white",
        xaxis_title="T",
        yaxis_title="W(L,T)",
        legend_title="Spatial size",
    )
    _write_figure_html(fig, path)


def analyze(
    input_paths: list[Path],
    output_dir: Path,
    thermalization_by_input: dict[Path, int],
    block_size: int,
    n_bootstrap: int,
    seed: int,
    *,
    run_dirs: list[Path] | None = None,
    group_signature: tuple[tuple[str, Any], ...] = tuple(),
    pair_coverage: dict[str, Any] | None = None,
) -> dict[str, Any]:
    values_by_loop, conf_range_by_loop = read_combined_gradient_flow_wtemp(
        input_paths,
        thermalization_by_input,
        label="final stats",
    )
    if not values_by_loop:
        raise ValueError("no Wilson-loop measurements remain after thermalization cut")

    output_dir.mkdir(parents=True, exist_ok=True)
    bootstrap_dir = output_dir / "bootstrap"
    bootstrap_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    bootstrap_plans: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray, int, int]] = {}

    records: list[dict[str, Any]] = []
    for flow_time, l_size, t_size in sorted(values_by_loop):
        values = np.asarray(values_by_loop[(flow_time, l_size, t_size)], dtype=float)
        plan = bootstrap_plans.get(int(values.size))
        if plan is None:
            plan = build_bootstrap_plan(int(values.size), block_size, n_bootstrap, rng)
            bootstrap_plans[int(values.size)] = plan
        samples, n_blocks, tail = bootstrap_means(values, plan)
        error = float(np.std(samples)) if samples.size > 0 else None

        flow_dir = bootstrap_dir / f"t_over_a2_{filename_token(flow_time)}"
        flow_dir.mkdir(parents=True, exist_ok=True)
        bootstrap_path = flow_dir / f"W_L_{l_size}_T_{t_size}.npy"
        np.save(bootstrap_path, samples)

        first_conf_id, last_conf_id = conf_range_by_loop[(flow_time, l_size, t_size)]
        records.append(
            {
                "t_over_a2": float(flow_time),
                "L": int(l_size),
                "T": int(t_size),
                "n_measurements": int(values.size),
                "first_conf_id": float(first_conf_id),
                "last_conf_id": float(last_conf_id),
                "block_size": int(block_size),
                "n_blocks": int(n_blocks),
                "tail": int(tail),
                "mean": float(np.mean(values, dtype=np.float64)),
                "error": error,
                "bootstrap_error": error,
                "bootstrap_path": relative_to_output(bootstrap_path, output_dir),
            }
        )

    summary = {
        "schema_version": 1,
        "input_paths": [str(input_path) for input_path in input_paths],
        "run_dirs": [str(run_dir) for run_dir in (run_dirs or [])],
        "group_signature": [
            {"name": name, "value": jsonable_value(value)}
            for name, value in group_signature
        ],
        "pair_coverage": pair_coverage,
        "thermalization_by_input": {
            str(input_path): int(thermalization_by_input[input_path])
            for input_path in input_paths
        },
        "block_size": int(block_size),
        "n_bootstrap": int(n_bootstrap),
        "seed": int(seed),
        "n_runs": len(run_dirs or input_paths),
        "n_wilson_loops": len(records),
        "flow_times": sorted({record["t_over_a2"] for record in records}),
        "records": records,
        "saved_at": datetime.now(timezone.utc).isoformat(),
    }

    summary_path = output_dir / "summary.json"
    with summary_path.open("w", encoding="utf-8") as stream:
        json.dump(summary, stream, indent=2, sort_keys=True)
        stream.write("\n")

    table_path = output_dir / "wilson_loop_stats.dat"
    with table_path.open("w", encoding="utf-8") as stream:
        stream.write(
            "# t_over_a2 L T n_measurements first_conf_id last_conf_id "
            "block_size n_blocks tail mean bootstrap_error bootstrap_path\n"
        )
        for record in records:
            stream.write(
                f"{record['t_over_a2']:.12g} {record['L']} {record['T']} "
                f"{record['n_measurements']} "
                f"{record['first_conf_id']:.12g} {record['last_conf_id']:.12g} "
                f"{record['block_size']} {record['n_blocks']} {record['tail']} "
                f"{record['mean']:.12g} {record['bootstrap_error']:.12g} "
                f"{record['bootstrap_path']}\n"
            )

    plot_path = output_dir / "wilson_loop_means.html"
    save_wilson_loop_plot(plot_path, records)

    return {
        "summary_path": str(summary_path),
        "table_path": str(table_path),
        "plot_path": str(plot_path),
        "bootstrap_dir": str(bootstrap_dir),
        "n_wilson_loops": len(records),
    }


def main() -> int:
    args = parse_args()
    if args.THERM is not None and args.THERM < 0:
        sys.exit("THERM must be non-negative")
    if args.BLOCK is not None and args.BLOCK < 1:
        sys.exit("BLOCK must be >= 1")
    if args.n_bootstrap < 1:
        sys.exit("--n-bootstrap must be >= 1")
    if args.min_group_size < 1:
        sys.exit("--min-group-size must be >= 1")

    try:
        date_min, date_max = resolve_date_filters(args)
        filter_criteria = parse_filter_tokens(args.filter)
        discovered, groups = resolve_input_groups(
            args.path,
            date_min,
            date_max,
            args.min_group_size,
            filter_criteria,
        )
    except OSError as exc:
        sys.exit(f"failed to read input: {exc}")
    except ValueError as exc:
        sys.exit(str(exc))

    if args.list_groups:
        diff_fields = differing_group_fields(groups)
        for index, group in enumerate(groups, start=1):
            print(f"[group {index}] {describe_group(group['signature'])}")
            diff_text = format_group_difference(group["signature"], diff_fields)
            if diff_text:
                print(f"  differs: {diff_text}")
            print_pair_coverage_if_mixed(group)
            for run_dir in group["run_dirs"]:
                print(f"  {run_dir}")
            if not group["run_dirs"]:
                for input_path in group["input_paths"]:
                    print(f"  {input_path}")
        return 0

    print("Selected analysis groups:")
    diff_fields = differing_group_fields(groups)
    for index, group in enumerate(groups, start=1):
        print(f"[group {index}] {describe_group(group['signature'])}")
        diff_text = format_group_difference(group["signature"], diff_fields)
        if diff_text:
            print(f"  differs: {diff_text}")
        print_pair_coverage_if_mixed(group)
        for run_dir in group["run_dirs"]:
            print(f"  {run_dir}")
        if not group["run_dirs"]:
            for input_path in group["input_paths"]:
                print(f"  {input_path}")

    base_input_path = normalize_path(args.path)
    if args.output_dir is not None:
        output_root = args.output_dir
    elif base_input_path.is_file():
        output_root = base_input_path.parent / "gradient_flow_wtemp_analysis"
    else:
        output_root = base_input_path / "gradient_flow_wtemp_analysis"

    results: list[dict[str, Any]] = []
    for group in groups:
        group_output_dir = output_root
        if discovered:
            group_output_dir = output_root / short_group_token(group["signature"])
        group_output_dir.mkdir(parents=True, exist_ok=True)

        if args.THERM is None:
            thermalization_plot_path = group_output_dir / "thermalization_preview.html"
            try:
                save_thermalization_preview_plot(
                    thermalization_plot_path,
                    group,
                    suggested_cut=int(run_evaluation.THERMALIZATION_STEPS),
                )
                print(f"Saved thermalization preview: {thermalization_plot_path}")
            except ValueError as exc:
                print(f"Could not save thermalization preview: {exc}")
            thermalization_by_input = prompt_thermalization_for_group(
                group,
                default_cut=int(run_evaluation.THERMALIZATION_STEPS),
            )
        else:
            thermalization_by_input = {
                input_path: int(args.THERM)
                for input_path in group["input_paths"]
            }
        save_decision_json(
            group_output_dir / "thermalization_choice.json",
            {
                "thermalization_by_input": {
                    str(input_path): int(thermalization_by_input[input_path])
                    for input_path in group["input_paths"]
                },
                "saved_at": datetime.now(timezone.utc).isoformat(),
            },
        )

        if args.BLOCK is None:
            default_block_size = int(run_evaluation.DEFAULT_BLOCK_SIZE)
            try:
                scan_values = block_size_scan_values(0, default_block_size)
                scan_records = build_bootstrap_block_size_scan_records(
                    group["input_paths"],
                    thermalization_by_input,
                    scan_values,
                    n_bootstrap=args.n_bootstrap,
                    seed=args.seed,
                )
                scan_json_path = group_output_dir / "bootstrap_block_size_scan.json"
                scan_plot_path = group_output_dir / "bootstrap_block_size_scan.html"
                save_decision_json(
                    scan_json_path,
                    {
                        "block_size_scan_values": scan_values,
                        "default_block_size": default_block_size,
                        "records": scan_records,
                        "saved_at": datetime.now(timezone.utc).isoformat(),
                    },
                )
                save_bootstrap_block_size_scan_plot(
                    scan_plot_path,
                    scan_records,
                    recommended_block_size=default_block_size,
                )
                print(f"Saved bootstrap block-size plot: {scan_plot_path}")
            except ValueError as exc:
                print(f"Could not save bootstrap block-size plot: {exc}")
            block_size = prompt_bootstrap_block_size(int(run_evaluation.DEFAULT_BLOCK_SIZE))
        else:
            block_size = int(args.BLOCK)
        save_decision_json(
            group_output_dir / "bootstrap_block_size_choice.json",
            {
                "block_size": int(block_size),
                "saved_at": datetime.now(timezone.utc).isoformat(),
            },
        )
        try:
            result = analyze(
                input_paths=group["input_paths"],
                output_dir=group_output_dir,
                thermalization_by_input=thermalization_by_input,
                block_size=block_size,
                n_bootstrap=args.n_bootstrap,
                seed=args.seed,
                run_dirs=group["run_dirs"],
                group_signature=group["signature"],
                pair_coverage=pair_coverage_summary(group["run_dirs"]) if group["run_dirs"] else None,
            )
        except OSError as exc:
            sys.exit(f"failed to write analysis output: {exc}")
        except ValueError as exc:
            sys.exit(str(exc))
        results.append(result)

    for result in results:
        print(f"Saved {result['n_wilson_loops']} Wilson-loop summaries")
        print(f"Summary: {result['summary_path']}")
        print(f"Table: {result['table_path']}")
        print(f"Plot: {result['plot_path']}")
        print(f"Bootstrap samples: {result['bootstrap_dir']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
