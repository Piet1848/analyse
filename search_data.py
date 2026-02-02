#!/usr/bin/env python3
"""
Search lattice runs by YAML parameters AND calculated results.
"""
from __future__ import annotations
import argparse
import os
import sys
from typing import Any, Dict, Tuple, List, Optional, get_type_hints, get_origin
import concurrent.futures

import run_evaluation
from load_input_yaml import load_params, MetropolisParams, GaugeObservableParams
import re

# --- Added tau_int and updated types ---
CALCULATED_FIELDS = {
    "W_R_T": float,
    "W_R_T_err": float,
    "V_R": dict,      # Changed to dict
    "V_R_err": dict,  # Changed to dict
    "r0": float,
    "r0_err": float,
    "a": float,
    "a_err": float,
    "tau_int": float, 
    "block_size": int,
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

def parse_bool(s: str) -> bool:
    sl = s.lower()
    if sl in ("true", "1", "yes", "y", "t"): return True
    if sl in ("false", "0", "no", "n", "f"): return False
    raise ValueError(f"Cannot parse boolean from '{s}'")

def parse_dynamic_token(tok: str):
    """Parses tokens like V_R5, W_R5_T3, W_R5_T, or W_R_T4. Handles _err suffix."""
    is_err = tok.endswith("_err")
    check_tok = tok[:-4] if is_err else tok

    # V_R<R>
    m = re.fullmatch(r"V_R(\d+)", check_tok)
    if m: return ("V_R_err" if is_err else "V_R"), {"R": int(m.group(1))}
    # W_R<R>_T<T>
    m = re.fullmatch(r"W_R(\d+)_T(\d+)", check_tok)
    if m and not is_err: return "W_R_T", {"R": int(m.group(1)), "T": int(m.group(2))}
    # W_R<R>_T
    m = re.fullmatch(r"W_R(\d+)_T", check_tok)
    if m and not is_err: return "W_R_T", {"R": int(m.group(1))}
    # W_R_T<T>
    m = re.fullmatch(r"W_R_T(\d+)", check_tok)
    if m and not is_err: return "W_R_T", {"T": int(m.group(1))}
    # creutz_P<R>
    m = re.fullmatch(r"creutz_P(\d+)", check_tok)
    if m: return ("creutz_P_err" if is_err else "creutz_P"), {"R": int(m.group(1))}
    # a_creutz<R>
    m = re.fullmatch(r"a_creutz(\d+)", check_tok)
    if m: return ("a_creutz_err" if is_err else "a_creutz"), {"R": int(m.group(1))}
    
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
        # For simple list/dict types, we treat input as string for now
        if typ is dict: return value_str 
        raise ValueError(f"List/tuple fields not supported: {typ}")
    
    if typ is int: return int(value_str)
    if typ is float: return float(value_str)
    if typ is bool: return parse_bool(value_str)
    if typ is dict: return value_str
    return value_str

def parse_tokens(tokens: list[str]) -> tuple[Dict[str, Any], list[str]]:
    criteria: Dict[str, Any] = {}
    outputs: list[str] = []
    for tok in tokens:
        if "=" in tok:
            name, value_str = tok.split("=", 1)
            name = name.strip(); value_str = value_str.strip()
            if name not in FIELD_MAP:
                raise KeyError(f"Unknown parameter '{name}'")
            _, typ = FIELD_MAP[name]
            if FIELD_MAP[name][0] == "calc":
                 print(f"Warning: Filtering by calculated field '{name}' triggers calculation.", file=sys.stderr)
            value = convert_value(value_str, typ)
            criteria[name] = value
        else:
            name = tok.strip()
            if not name: continue
            info = get_field_info(name)
            if not info:
                raise KeyError(f"Unknown parameter '{name}'")
            outputs.append(name)
    return criteria, outputs

def matches_criteria(metro: MetropolisParams, gauge: GaugeObservableParams, criteria: Dict[str, Any]) -> bool:
    for name, expected in criteria.items():
        block_name, _ = FIELD_MAP[name]
        if block_name == "calc": continue 
        obj = metro if block_name == "metro" else gauge
        actual = getattr(obj, name)
        if actual != expected: return False
    return True

def find_matching_runs(root: str, criteria: Dict[str, Any]) -> list[str]:
    matches: list[str] = []
    for dirpath, dirnames, filenames in os.walk(root):
        if "input.yaml" not in filenames: continue
        yaml_path = os.path.join(dirpath, "input.yaml")
        try:
            metro, gauge = load_params(yaml_path)
        except Exception: continue
        if matches_criteria(metro, gauge, criteria):
            matches.append(dirpath)
    return matches

def process_row(path: str, outputs: List[str]) -> Tuple[str, List[Any]]:
    try:
        metro, gauge = load_params(os.path.join(path, "input.yaml"))
        needs_calc = any(get_field_info(out)[0] == "calc" for out in outputs)
        calc_data = run_evaluation.get_or_calculate(path) if needs_calc else None

        vals = []
        for name in outputs:
            info = get_field_info(name)
            block, _ = info
            if block == "calc":
                if not calc_data or "error" in calc_data:
                    vals.append(None)
                elif name in calc_data:
                    vals.append(calc_data[name])
                else:
                    # Dynamic extraction
                    base, params = parse_dynamic_token(name)
                    if base in ["V_R", "V_R_err", "creutz_P", "creutz_P_err", "a_creutz", "a_creutz_err"]:
                        vals.append(calc_data.get(base, {}).get(str(params["R"])))
                    elif base.startswith("W_R_T"):
                        w_dict = calc_data.get("W_R_T", {})
                        r, t = params.get("R"), params.get("T")
                        if r is not None and t is not None:
                            vals.append(w_dict.get(f"{r},{t}"))
                        elif r is not None:
                            # Return dict of all T for this R
                            sub = {k.split(',')[1]: v for k, v in w_dict.items() if k.startswith(f"{r},")}
                            vals.append(sub if sub else None)
                        elif t is not None:
                            # Return dict of all R for this T
                            sub = {k.split(',')[0]: v for k, v in w_dict.items() if k.endswith(f",{t}")}
                            vals.append(sub if sub else None)
                        else: vals.append(None)
            else:
                obj = metro if block == "metro" else gauge
                vals.append(getattr(obj, name))
        return path, vals
    except Exception:
        return path, [None] * len(outputs)

def main() -> None:
    parser = argparse.ArgumentParser(description="Search data runs and calculate missing values.")
    parser.add_argument("root", help="Root data directory")
    parser.add_argument("tokens", nargs="*", help="Filters (NAME=VAL) and outputs (NAME)")
    parser.add_argument("--workers", type=int, default=2, help="Number of parallel workers")
    args = parser.parse_args()

    if not args.tokens:
        print("Available Parameters:")
        for k, v in FIELD_MAP.items():
            if v[0] != "calc": print(f"  {k}")
        print("\nCalculated Fields:")
        for k, v in FIELD_MAP.items():
            if v[0] == "calc": print(f"  {k}")
        return

    try:
        criteria, outputs = parse_tokens(args.tokens)
    except Exception as e:
        sys.exit(f"Error: {e}")

    matches = find_matching_runs(args.root, criteria)
    if not matches:
        print("No matching runs found.")
        return

    print(f"Found {len(matches)} matching runs. Processing with {args.workers} workers...", file=sys.stderr)

    rows = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(process_row, path, outputs): path for path in matches}
        for i, future in enumerate(concurrent.futures.as_completed(futures)):
            path, vals = future.result()
            rows.append((path, vals))
            print(f"Processed {i+1}/{len(matches)}", file=sys.stderr, end="\r")

    print(" " * 60, file=sys.stderr, end="\r")
    
    # Sort
    def sort_key(row):
        key_parts = []
        for val in row[1]:
            if val is None:
                key_parts.append((0, ""))
            elif isinstance(val, (dict, list)):
                # Convert dicts/lists to string for sorting purposes
                key_parts.append((1, str(val)))
            else:
                key_parts.append((1, val))
        return tuple(key_parts)
        
    rows.sort(key=sort_key)

    print("\t".join(["path"] + outputs))
    for path, vals in rows:
        out_strs = [str(v) if v is not None else "N/A" for v in vals]
        print("\t".join([path] + out_strs))

if __name__ == "__main__":
    main()