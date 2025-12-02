#!/usr/bin/env python3
"""
Search lattice runs by YAML parameters.

Usage examples
--------------
# From ~/masterarbeit (where ./data exists):
python3 search_data.py ./data epsilon1=-0.75 beta epsilon1

# Explicit root and single condition:
python3 search_data.py ./data nSweep=10000 beta
"""
from __future__ import annotations

import argparse
import os
import sys
from typing import Any, Dict, Tuple, List, Optional, get_type_hints, get_origin

# Adjust this import if your file has a different name
from load_input_yaml import load_params, MetropolisParams, GaugeObservableParams


def build_field_map() -> Dict[str, Tuple[str, Any]]:
    """
    Map parameter name -> (block_name, type).

    block_name is "metro" or "gauge".
    """
    field_map: Dict[str, Tuple[str, Any]] = {}

    metro_types = get_type_hints(MetropolisParams)
    for name, typ in metro_types.items():
        field_map[name] = ("metro", typ)

    gauge_types = get_type_hints(GaugeObservableParams)
    for name, typ in gauge_types.items():
        field_map[name] = ("gauge", typ)

    return field_map


FIELD_MAP = build_field_map()


def parse_bool(s: str) -> bool:
    sl = s.lower()
    if sl in ("true", "1", "yes", "y", "t"):
        return True
    if sl in ("false", "0", "no", "n", "f"):
        return False
    raise ValueError(f"Cannot parse boolean from '{s}'")


def convert_value(value_str: str, typ: Any) -> Any:
    """
    Convert command line string to the appropriate type based on the
    dataclass annotation.
    """
    origin = get_origin(typ)

    # Only support scalar fields for now
    if origin is not None:
        raise ValueError(f"List/tuple fields not supported in search: type={typ}")

    if typ is int:
        return int(value_str)
    if typ is float:
        return float(value_str)
    if typ is bool:
        return parse_bool(value_str)

    # Fallback: compare as string
    return value_str


def parse_tokens(tokens: list[str]) -> tuple[Dict[str, Any], list[str]]:
    """
    Split tokens into:
     - criteria dict from tokens containing '=' (NAME=VALUE)
     - outputs list from tokens without '=' (names to print for each match)

    Returns (criteria, outputs)
    """
    criteria: Dict[str, Any] = {}
    outputs: list[str] = []

    for tok in tokens:
        if "=" in tok:
            name, value_str = tok.split("=", 1)
            name = name.strip()
            value_str = value_str.strip()

            if name not in FIELD_MAP:
                available = ", ".join(sorted(FIELD_MAP.keys()))
                raise KeyError(
                    f"Unknown parameter '{name}'. Known parameters are: {available}"
                )

            _, typ = FIELD_MAP[name]
            value = convert_value(value_str, typ)
            criteria[name] = value
        else:
            name = tok.strip()
            if name == "":
                continue
            if name not in FIELD_MAP:
                available = ", ".join(sorted(FIELD_MAP.keys()))
                raise KeyError(
                    f"Unknown output parameter '{name}'. Known parameters are: {available}"
                )
            outputs.append(name)

    return criteria, outputs


def matches_criteria(metro: MetropolisParams,
                     gauge: GaugeObservableParams,
                     criteria: Dict[str, Any]) -> bool:
    """
    Check whether given params satisfy all criteria.
    """
    for name, expected in criteria.items():
        block_name, _ = FIELD_MAP[name]
        obj = metro if block_name == "metro" else gauge
        actual = getattr(obj, name)

        if actual != expected:
            return False
    return True


def find_matching_runs(root: str, criteria: Dict[str, Any]) -> list[str]:
    """
    Walk 'root', look for input.yaml in subfolders, apply criteria,
    and return a list of directories with matches.
    """
    matches: list[str] = []

    for dirpath, dirnames, filenames in os.walk(root):
        if "input.yaml" not in filenames:
            continue

        yaml_path = os.path.join(dirpath, "input.yaml")
        try:
            metro, gauge = load_params(yaml_path)
        except Exception as e:
            print(f"Warning: failed to load '{yaml_path}': {e}", file=sys.stderr)
            continue

        if matches_criteria(metro, gauge, criteria):
            matches.append(dirpath)

    return matches


def get_value_str(name: str, metro: MetropolisParams, gauge: GaugeObservableParams) -> str:
    """
    Return a string representation for the given parameter name taken from
    metro or gauge dataclass.
    """
    block_name, _ = FIELD_MAP[name]
    obj = metro if block_name == "metro" else gauge
    val = getattr(obj, name)
    # Basic formatting for booleans and None; other types use str()
    if isinstance(val, bool):
        return "true" if val else "false"
    if val is None:
        return ""
    return str(val)


def print_available_parameters():
    print("Available parameters:\n")

    print("MetropolisParams:")
    for name, (block, _) in FIELD_MAP.items():
        if block == "metro":
            print(f"  {name}")
    print()

    print("GaugeObservableParams:")
    for name, (block, _) in FIELD_MAP.items():
        if block == "gauge":
            print(f"  {name}")
    print()

    from typing import List, Optional

def search_paths(root: str,
                 filters: Optional[Dict[str, Any]] = None,
                 filter_tokens: Optional[List[str]] = None) -> List[str]:
    """
    Programmatic API returning a list of matching run directories.

    - root: path to data root (same as CLI first argument)
    - filters: optional dict like {'epsilon1': -0.75, 'beta': 2.3}
    - filter_tokens: optional list of CLI-style filter tokens like ['epsilon1=-0.75','beta=2.3'].
                     If provided, it is parsed with parse_tokens() and overrides `filters`.

    Returns: list of matching directory paths (unsorted; same order as os.walk).
    """
    if filter_tokens is not None:
        parsed_filters, _ = parse_tokens(filter_tokens)
        filters = parsed_filters

    if filters is None:
        filters = {}

    # validate filter keys
    for k in filters.keys():
        if k not in FIELD_MAP:
            raise KeyError(f"Unknown filter parameter '{k}'")

    return find_matching_runs(root, filters)



def main() -> None:
    parser = argparse.ArgumentParser(
        description="Search data runs by parameters in input.yaml."
    )
    parser.add_argument(
        "root",
        help="Root data directory (e.g. ./data or ./data/20251124).",
    )
    parser.add_argument(
        "tokens",
        nargs="*",
        help="Filters and outputs. Filters must be NAME=VALUE. Non-equal tokens are output field names.",
    )

    args = parser.parse_args()

    # If no tokens → print full list and exit
    if len(args.tokens) == 0:
        print_available_parameters()
        return

    # Split tokens into filters (NAME=VALUE) and outputs (names)
    try:
        criteria, outputs = parse_tokens(args.tokens)
    except (ValueError, KeyError) as e:
        print(f"Error parsing arguments: {e}", file=sys.stderr)
        sys.exit(1)

    matches = find_matching_runs(args.root, criteria)

    if not matches:
        print("No matching runs found.")
        return

    # Read raw values for outputs (typed) for a given path
    def read_output_values_raw(path):
        metro, gauge = load_params(os.path.join(path, "input.yaml"))
        vals = []
        for name in outputs:
            block, _ = FIELD_MAP[name]
            obj = metro if block == "metro" else gauge
            vals.append(getattr(obj, name))
        return metro, gauge, vals

    rows = []
    for path in matches:
        if outputs:
            try:
                metro, gauge, vals = read_output_values_raw(path)
            except Exception as e:
                print(f"Warning: failed to load '{path}/input.yaml': {e}", file=sys.stderr)
                continue
            rows.append((path, metro, gauge, vals))
        else:
            rows.append((path, None, None, []))

    # Sorting
    if outputs:
        type_list = [FIELD_MAP[name][1] for name in outputs]

        def convert_for_sort(v, typ):
            # handle None
            if v is None:
                return (0, "")  # ensure None sorts before real values
            # basic scalar typing as in convert_value
            origin = get_origin(typ)
            if origin is not None:
                # fallback to string for container/complex types
                return (1, str(v))
            try:
                if typ is int:
                    return (1, int(v))
                if typ is float:
                    return (1, float(v))
                if typ is bool:
                    return (1, bool(v))
            except Exception:
                # conversion failed, fall back to string
                return (1, str(v))
            return (1, str(v))

        def sort_key(entry):
            _, _, _, vals = entry
            return tuple(convert_for_sort(val, typ) for val, typ in zip(vals, type_list))

        rows.sort(key=sort_key)
    else:
        rows.sort(key=lambda r: r[0])

    # Print header (tab separated)
    header = ["path"] + outputs
    print("\t".join(header))

    # Print rows: format values using get_value_str for consistent formatting
    for path, metro, gauge, vals in rows:
        if outputs:
            out_strs = []
            for name, raw in zip(outputs, vals):
                # use get_value_str which loads from metro/gauge; to avoid reloading,
                # create a small helper that formats raw values similarly:
                if isinstance(raw, bool):
                    out_strs.append("true" if raw else "false")
                elif raw is None:
                    out_strs.append("")
                else:
                    out_strs.append(str(raw))
            print("\t".join([path] + out_strs))
        else:
            print(path)


if __name__ == "__main__":
    main()
