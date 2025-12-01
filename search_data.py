#!/usr/bin/env python3
"""
Search lattice runs by YAML parameters.

Usage examples
--------------
# From ~/masterarbeit (where ./data exists):
python find_runs.py ./data beta=2.3 epsilon1=-0.75

# Explicit root and single condition:
python find_runs.py ./data nSweep=10000
"""

import argparse
import os
import sys
from typing import Any, Dict, Tuple, get_type_hints, get_origin

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


def parse_criteria(raw_criteria: list[str]) -> Dict[str, Any]:
    """
    Parse a list like ["beta=2.3", "epsilon1=-0.75"] into a dict
    { "beta": 2.3, "epsilon1": -0.75 } with proper types.
    """
    criteria: Dict[str, Any] = {}

    for item in raw_criteria:
        if "=" not in item:
            raise ValueError(f"Invalid criterion '{item}', expected NAME=VALUE")

        name, value_str = item.split("=", 1)
        name = name.strip()
        value_str = value_str.strip()

        if name not in FIELD_MAP:
            available = ", ".join(sorted(FIELD_MAP.keys()))
            raise KeyError(
                f"Unknown parameter '{name}'. "
                f"Known parameters are: {available}"
            )

        _, typ = FIELD_MAP[name]
        value = convert_value(value_str, typ)
        criteria[name] = value

    return criteria


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


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Search data runs by parameters in input.yaml."
    )
    parser.add_argument(
        "root",
        help="Root data directory (e.g. ./data or ./data/20251124).",
    )
    parser.add_argument(
        "criteria",
        nargs="*",
        help="Parameter filters of the form NAME=VALUE, e.g. beta=2.3 epsilon1=-0.75",
    )

    args = parser.parse_args()
    
    # If no parameters → print full list and exit
    if len(args.criteria) == 0:
        print_available_parameters()
        return

    try:
        criteria = parse_criteria(args.criteria)
    except (ValueError, KeyError) as e:
        print(f"Error parsing criteria: {e}", file=sys.stderr)
        sys.exit(1)

    matches = find_matching_runs(args.root, criteria)

    if not matches:
        print("No matching runs found.")
        return

    print("Matching run directories:")
    for path in sorted(matches):
        print(path)


if __name__ == "__main__":
    main()
