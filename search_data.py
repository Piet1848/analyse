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

    return None, None


def get_field_info(name: str):
    """Unified lookup for static and dynamic fields."""
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
            name = name.strip()
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
    if name in calc_data:
        return calc_data[name]

    base, params = parse_dynamic_token(name)
    if not base:
        return None

    if base in ["V_R", "V_R_err", "creutz_P", "creutz_P_err", "a_creutz", "a_creutz_err"]:
        return calc_data.get(base, {}).get(str(params["R"]))

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
                vals.append(getattr(obj, name))
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
) -> Tuple[str, Dict[str, Any]]:
    run_evaluation = _get_run_evaluation()
    return calc_key, run_evaluation.get_or_calculate(
        path,
        force_recalc=force_calculate,
        calc_workers=calc_workers,
        load_workers=load_workers,
        combine_equivalent_runs=combine_equivalent_runs,
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

    print("	".join(["path"] + outputs))
    for row_label, vals in rows:
        out_strs = [str(v) if v is not None else "N/A" for v in vals]
        print("	".join([row_label] + out_strs))


if __name__ == "__main__":
    main()
