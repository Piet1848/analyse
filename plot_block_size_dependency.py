#!/usr/bin/env python3
"""
Plot the bootstrap error of a scalar Calculator variable as a function of
bootstrap blocking size.

The script intentionally bypasses cached analysis output and recalculates the
requested variable for every block size. Data loading and run grouping follow
the same logic as run_evaluation.py.

Heavy numeric kernels are delegated to NumPy/BLAS/OpenMP threading when
available. The script itself does not create its own worker threads.
"""
from __future__ import annotations

import argparse
import csv
import inspect
import math
import os
import re
import sys
from pathlib import Path
from typing import Any

NUMPY_THREAD_ENV_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "BLIS_NUM_THREADS",
)


def apply_numpy_thread_settings_from_argv(argv: list[str]) -> int | None:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--numpy-threads", type=int)
    args, _ = parser.parse_known_args(argv[1:])
    if args.numpy_threads is None:
        return None

    thread_count = max(1, int(args.numpy_threads))
    for env_name in NUMPY_THREAD_ENV_VARS:
        os.environ[env_name] = str(thread_count)
    return thread_count


CONFIGURED_NUMPY_THREADS = apply_numpy_thread_settings_from_argv(sys.argv)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from finalized_analysis_helpers import save_matplotlib_figure
import run_evaluation
from calculator import Calculator
from load_input_yaml import load_params


def parse_bool(text: str) -> bool:
    lowered = text.strip().lower()
    if lowered in {"1", "true", "t", "yes", "y", "on"}:
        return True
    if lowered in {"0", "false", "f", "no", "n", "off"}:
        return False
    raise ValueError(f"Cannot parse boolean from {text!r}")


def parse_value(text: str) -> Any:
    stripped = text.strip()
    lowered = stripped.lower()

    if lowered in {"none", "null"}:
        return None

    try:
        return parse_bool(stripped)
    except ValueError:
        pass

    try:
        return int(stripped)
    except ValueError:
        pass

    try:
        return float(stripped)
    except ValueError:
        pass

    return stripped


def parse_param_tokens(tokens: list[str]) -> dict[str, Any]:
    params: dict[str, Any] = {}
    for token in tokens:
        if "=" not in token:
            raise ValueError(f"Expected NAME=VALUE token, got {token!r}")
        name, raw_value = token.split("=", 1)
        name = name.strip()
        if not name:
            raise ValueError(f"Invalid empty parameter name in token {token!r}")
        params[name] = parse_value(raw_value)
    return params


def parse_block_sizes(args: argparse.Namespace) -> list[int]:
    if args.block_sizes:
        values: list[int] = []
        for chunk in args.block_sizes.split(","):
            stripped = chunk.strip()
            if not stripped:
                continue
            values.append(int(stripped))
    else:
        if args.block_step < 1:
            raise ValueError("--block-step must be at least 1")
        values = list(range(args.min_block_size, args.max_block_size + 1, args.block_step))

    block_sizes = sorted({int(v) for v in values if int(v) >= 1})
    if not block_sizes:
        raise ValueError("No valid block sizes selected.")
    return block_sizes


def ensure_variable_exists(variable: str) -> None:
    if not Calculator._registry:
        Calculator._populate_registry()
    if variable not in Calculator._registry:
        available = ", ".join(sorted(Calculator._registry))
        raise KeyError(f"Unknown variable {variable!r}. Available variables: {available}")


def default_variable_params(variable: str, input_dir: Path) -> dict[str, Any]:
    params: dict[str, Any] = {}

    if variable == "V_R":
        params["t_min"] = run_evaluation.DEFAULT_V_R_T_MIN
        params["t_max"] = run_evaluation.DEFAULT_V_R_T_MAX
    elif variable == "r0":
        params.update(
            {
                "t_min": run_evaluation.DEFAULT_R0_T_MIN,
                "t_max": run_evaluation.DEFAULT_R0_T_MAX,
                "target_force": run_evaluation.DEFAULT_SOMMER_TARGET,
                "r_min": run_evaluation.DEFAULT_R0_R_MIN,
            }
        )
    elif variable == "F_chi":
        params["t_large"] = run_evaluation.DEFAULT_R0_CHI_T_LARGE
    elif variable == "r0_chi":
        params.update(
            {
                "t_large": run_evaluation.DEFAULT_R0_CHI_T_LARGE,
                "target_force": run_evaluation.DEFAULT_SOMMER_TARGET,
                "max_rel_err": run_evaluation.DEFAULT_R0_CHI_MAX_REL_ERR,
                "use_weighted_fit": run_evaluation.DEFAULT_R0_CHI_USE_WEIGHTED_FIT,
                "fit_window": run_evaluation.DEFAULT_R0_CHI_FIT_WINDOW,
                "discard_negative": run_evaluation.DEFAULT_R0_CHI_DISCARD_NEGATIVE,
                "r_min": run_evaluation.DEFAULT_R0_CHI_R_MIN,
            }
        )
    elif variable == "volume_r0":
        metro, _ = load_params(str(input_dir / "input.yaml"))
        params.update(
            {
                "L0": metro.L0,
                "L1": metro.L1,
                "L2": metro.L2,
                "L3": metro.L3,
                "r0_phys": 0.5,
                "t_min": run_evaluation.DEFAULT_R0_T_MIN,
                "t_max": run_evaluation.DEFAULT_R0_T_MAX,
                "target_force": run_evaluation.DEFAULT_SOMMER_TARGET,
                "r_min": run_evaluation.DEFAULT_R0_R_MIN,
            }
        )

    return params


def build_variable_params(variable: str, user_params: dict[str, Any], input_dir: Path) -> dict[str, Any]:
    params = default_variable_params(variable, input_dir)
    params.update(user_params)
    return params


def load_file_data(
    path: str,
    combine_equivalent_runs: bool,
    load_workers: int,
    verbose: bool,
):
    abs_path = os.path.abspath(path)
    if not os.path.isdir(abs_path):
        raise FileNotFoundError(f"Directory not found: {abs_path}")

    if combine_equivalent_runs:
        analysis_id, grouped_paths = run_evaluation._discover_equivalent_runs(abs_path)
    else:
        analysis_id = f"run_{run_evaluation.get_run_id(abs_path)}"
        grouped_paths = [abs_path]

    prefix = f"[{os.path.basename(abs_path)}]"
    file_data, aggregation = run_evaluation._load_combined_w_temp(
        grouped_paths,
        verbose=verbose,
        prefix=prefix,
        load_workers=load_workers,
    )
    if file_data is None:
        raise RuntimeError("No W_temp data found for the selected run(s).")

    return abs_path, analysis_id, grouped_paths, aggregation, file_data


def compute_tau_hint(file_data) -> tuple[float | None, int | None]:
    if not any(obs.name in {"plaquette", "retrace"} for obs in file_data.observables):
        return None, None

    calc = Calculator(file_data, n_bootstrap=run_evaluation.DEFAULT_N_BOOTSTRAP)
    try:
        tau = calc.get_variable("tau_int", obs_name="plaquette").get()
    except KeyError:
        try:
            tau = calc.get_variable("tau_int", obs_name="retrace").get()
        except KeyError:
            return None, None

    if tau is None or not np.isfinite(tau):
        return None, None
    return float(tau), max(1, int(math.ceil(2 * float(tau))))


def to_float(value: Any) -> float:
    if value is None:
        return float("nan")
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"Expected a scalar numeric result, got {value!r}") from exc


def infer_required_w_rt_pairs(
    calc: Calculator,
    variable: str,
    params: dict[str, Any],
) -> list[tuple[int, int]]:
    available_pairs = set(calc.get_available_pairs())

    def unique_existing(pairs: list[tuple[int, int]]) -> list[tuple[int, int]]:
        seen: set[tuple[int, int]] = set()
        result: list[tuple[int, int]] = []
        for r_val, t_val in pairs:
            pair = (int(r_val), int(t_val))
            if pair in seen or pair not in available_pairs:
                continue
            seen.add(pair)
            result.append(pair)
        return result

    if variable == "W_R_T":
        return unique_existing([(int(params["R"]), int(params["T"]))])

    if variable == "V_R":
        r_val = int(params["R"])
        t_min = int(params.get("t_min", 1))
        t_max = params.get("t_max")
        return unique_existing(
            [(r_val, t_val) for t_val in calc.get_unique_Ts(R=r_val) if t_val >= t_min and (t_max is None or t_val <= int(t_max))]
        )

    if variable == "effective_mass":
        r_val = int(params["R"])
        t_val = int(params.get("T", 1))
        return unique_existing([(r_val, t_val), (r_val, t_val + 1)])

    if variable == "chi":
        r_int = int(float(params["R"]) - 0.5)
        t_int = int(float(params["T"]) - 0.5)
        return unique_existing(
            [
                (r_int, t_int),
                (r_int, t_int + 1),
                (r_int + 1, t_int),
                (r_int + 1, t_int + 1),
            ]
        )

    if variable == "F_chi":
        r_int = int(float(params["R"]) - 0.5)
        t_large = int(params["t_large"])
        return unique_existing(
            [
                (r_int, t_large),
                (r_int, t_large + 1),
                (r_int + 1, t_large),
                (r_int + 1, t_large + 1),
            ]
        )

    if variable in {"creutz_P", "a_creutz"}:
        r_val = int(params["R"])
        r_half = r_val // 2
        return unique_existing([(r_val, r_val), (r_half, r_half), (r_val, r_half), (r_half, r_val)])

    if variable in {"r0", "volume_r0"}:
        t_min = int(params.get("t_min", run_evaluation.DEFAULT_R0_T_MIN))
        t_max = params.get("t_max")
        r_min = int(params.get("r_min", run_evaluation.DEFAULT_R0_R_MIN))
        pairs = [
            (r_val, t_val)
            for r_val in calc.get_unique_Rs(r_min=r_min)
            for t_val in calc.get_unique_Ts(R=r_val)
            if t_val >= t_min and (t_max is None or t_val <= int(t_max))
        ]
        return unique_existing(pairs)

    if variable == "r0_chi":
        t_large = int(params.get("t_large", run_evaluation.DEFAULT_R0_CHI_T_LARGE))
        r_min = int(params.get("r_min", run_evaluation.DEFAULT_R0_CHI_R_MIN))
        unique_ls = calc.get_unique_Rs(r_min=r_min)
        pairs: list[tuple[int, int]] = []
        for l_val in unique_ls:
            if (l_val + 1) not in {r for r, _ in available_pairs}:
                continue
            pairs.extend(
                [
                    (l_val, t_large),
                    (l_val, t_large + 1),
                    (l_val + 1, t_large),
                    (l_val + 1, t_large + 1),
                ]
            )
        return unique_existing(pairs)

    return []


def calculate_one_block_size(
    file_data,
    variable: str,
    params: dict[str, Any],
    block_size: int,
    n_bootstrap: int,
) -> tuple[float, float, str]:
    try:
        calc = Calculator(file_data, n_bootstrap=n_bootstrap, step_size=block_size)
        calc.prime_w_rt_cache(pairs=infer_required_w_rt_pairs(calc, variable, params))
        variable_data = calc.get_variable(variable, **params)
        value = to_float(variable_data.get())
        estimate_error = to_float(variable_data.err()) if variable_data.err() is not None else float("nan")
        return value, estimate_error, "ok"
    except Exception as exc:
        return float("nan"), float("nan"), str(exc)


def calculate_series(
    file_data,
    variable: str,
    params: dict[str, Any],
    block_sizes: list[int],
    n_bootstrap: int,
) -> dict[str, np.ndarray]:
    block_size_array = np.asarray(block_sizes, dtype=int)
    values = np.full(block_size_array.shape, np.nan, dtype=float)
    estimate_errors = np.full(block_size_array.shape, np.nan, dtype=float)
    statuses = np.empty(block_size_array.shape, dtype=object)

    for idx, block_size in enumerate(block_size_array):
        value, estimate_error, status = calculate_one_block_size(
            file_data,
            variable,
            params,
            int(block_size),
            n_bootstrap,
        )
        values[idx] = value
        estimate_errors[idx] = estimate_error
        statuses[idx] = status

    return {
        "block_size": block_size_array,
        "value": values,
        "estimate_error": estimate_errors,
        "status": statuses,
    }


def sanitize_filename(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("_") or "plot"


def format_variable_label(variable: str, params: dict[str, Any]) -> str:
    if not params:
        return variable
    param_text = ", ".join(f"{key}={value}" for key, value in sorted(params.items()))
    return f"{variable}({param_text})"


def iter_rows(results: dict[str, np.ndarray]):
    for idx in range(len(results["block_size"])):
        yield {
            "block_size": int(results["block_size"][idx]),
            "value": float(results["value"][idx]),
            "estimate_error": float(results["estimate_error"][idx]),
            "status": str(results["status"][idx]),
        }


def write_csv(path: Path, results: dict[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["block_size", "value", "estimate_error", "status"],
        )
        writer.writeheader()
        writer.writerows(iter_rows(results))


def write_plot(
    path: Path,
    results: dict[str, np.ndarray],
    variable_label: str,
    run_label: str,
    recommended_block_size: int | None,
) -> dict[str, Path]:
    valid_mask = np.isfinite(results["estimate_error"])
    if not np.any(valid_mask):
        raise RuntimeError("All calculations failed; no finite error estimates available to plot.")

    x = results["block_size"][valid_mask].astype(float, copy=False)
    y = results["estimate_error"][valid_mask]

    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    ax.plot(x, y, "o", markersize=1.5, alpha=0.6)

    if recommended_block_size is not None:
        ax.axvline(
            recommended_block_size,
            color="0.4",
            linestyle="--",
            linewidth=1.0,
            label=f"tau-based block size = {recommended_block_size}",
        )
        ax.legend()

    ax.set_xlabel("Blocking size")
    ax.set_ylabel(f"Bootstrap error of {variable_label}")
    ax.set_title(f"Bootstrap error vs blocking size\n{run_label}")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    written = save_matplotlib_figure(fig, path, dpi=160)
    plt.close(fig)
    if not written:
        raise RuntimeError(f"No plot exports could be written for {path}.")
    return written


def print_table(results: dict[str, np.ndarray]) -> None:
    print("block_size\tvalue\testimate_error\tstatus")
    for idx in range(len(results["block_size"])):
        value_float = float(results["value"][idx])
        estimate_error_float = float(results["estimate_error"][idx])
        value = "nan" if not np.isfinite(value_float) else f"{value_float:.10g}"
        estimate_error = "nan" if not np.isfinite(estimate_error_float) else f"{estimate_error_float:.10g}"
        print(
            f"{int(results['block_size'][idx])}\t{value}\t{estimate_error}\t{results['status'][idx]}"
        )


def build_output_path(
    explicit_path: str | None,
    analysis_id: str,
    variable: str,
    suffix: str,
) -> Path:
    if explicit_path:
        return Path(explicit_path)
    base_name = f"{analysis_id}__{sanitize_filename(variable)}__block_size"
    return Path(f"{base_name}{suffix}")


def describe_signature(variable: str) -> str:
    func = Calculator._registry[variable]
    signature = inspect.signature(func)
    params = []
    for name, param in signature.parameters.items():
        if name == "self":
            continue
        params.append(str(param))
    return ", ".join(params)


def numpy_thread_summary(args: argparse.Namespace) -> str:
    if args.numpy_threads is not None:
        return str(max(1, int(args.numpy_threads)))
    if CONFIGURED_NUMPY_THREADS is not None:
        return str(CONFIGURED_NUMPY_THREADS)
    for env_name in NUMPY_THREAD_ENV_VARS:
        value = os.environ.get(env_name)
        if value:
            return f"env:{env_name}={value}"
    return "library default"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Plot the bootstrap error of a Calculator variable as a function of blocking size.",
    )
    parser.add_argument("path", help="Run directory containing input.yaml and W_temp.out")
    parser.add_argument("variable", help="Calculator variable name, for example V_R or r0")
    parser.add_argument(
        "params",
        nargs="*",
        help="Variable parameters as NAME=VALUE tokens, for example R=4 T=3",
    )
    parser.add_argument(
        "--block-sizes",
        help="Comma-separated block sizes. Overrides --min/--max.",
    )
    parser.add_argument("--min", "--min-block-size", dest="min_block_size", type=int, default=1, help="Minimum block size when using a range.")
    parser.add_argument("--max", "--max-block-size", dest="max_block_size", type=int, default=16, help="Maximum block size when using a range.")
    parser.add_argument("--block-step", type=int, default=1, help="Step for generated block sizes.")
    parser.add_argument(
        "--n-bootstrap",
        type=int,
        default=run_evaluation.DEFAULT_N_BOOTSTRAP,
        help="Number of bootstrap replicas per block size.",
    )
    parser.add_argument(
        "--numpy-threads",
        type=int,
        default=CONFIGURED_NUMPY_THREADS,
        help="Set NumPy/BLAS/OpenMP thread count for heavy numeric kernels before imports.",
    )
    parser.add_argument(
        "--load-workers",
        type=int,
        default=1,
        help="Files to load in parallel while combining equivalent runs.",
    )
    parser.add_argument(
        "--no-combine-equivalent-runs",
        action="store_true",
        help="Use only the selected run instead of combining equivalent runs.",
    )
    parser.add_argument("--output", help="Output image path. Default: <analysis_id>__<variable>__block_size.png")
    parser.add_argument("--csv", help="Output CSV path. Default: <analysis_id>__<variable>__block_size.csv")
    parser.add_argument("--verbose", action="store_true", help="Print data-loading progress to stderr.")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    try:
        ensure_variable_exists(args.variable)
        user_params = parse_param_tokens(args.params)
        block_sizes = parse_block_sizes(args)
    except Exception as exc:
        sys.exit(f"Error: {exc}")

    try:
        abs_path, analysis_id, grouped_paths, _, file_data = load_file_data(
            args.path,
            combine_equivalent_runs=not args.no_combine_equivalent_runs,
            load_workers=args.load_workers,
            verbose=args.verbose,
        )
    except Exception as exc:
        sys.exit(f"Error: {exc}")

    try:
        variable_params = build_variable_params(args.variable, user_params, Path(abs_path))
    except Exception as exc:
        sys.exit(f"Error while preparing variable parameters: {exc}")

    output_path = build_output_path(args.output, analysis_id, args.variable, ".png")
    csv_path = build_output_path(args.csv, analysis_id, args.variable, ".csv")
    variable_label = format_variable_label(args.variable, variable_params)

    rows = calculate_series(
        file_data=file_data,
        variable=args.variable,
        params=variable_params,
        block_sizes=block_sizes,
        n_bootstrap=args.n_bootstrap,
    )

    try:
        tau_hint, recommended_block_size = compute_tau_hint(file_data)
    except Exception:
        tau_hint, recommended_block_size = None, None

    write_csv(csv_path, rows)
    try:
        plot_paths = write_plot(
            output_path,
            rows,
            variable_label=variable_label,
            run_label=analysis_id,
            recommended_block_size=recommended_block_size,
        )
    except Exception as exc:
        print_table(rows)
        sys.exit(f"Error while writing plot: {exc}")

    print_table(rows)
    print()
    print(f"Loaded run path: {abs_path}")
    print(f"Analysis id: {analysis_id}")
    print(f"Combined runs: {len(grouped_paths)}")
    print(f"Variable signature: {args.variable}({describe_signature(args.variable)})")
    print(f"NumPy thread setting: {numpy_thread_summary(args)}")
    if tau_hint is None:
        print("Tau hint: unavailable from loaded data")
    else:
        print(f"Tau hint: tau_int={tau_hint:.6g}, recommended block size={recommended_block_size}")
    for fmt, path in sorted(plot_paths.items()):
        print(f"Plot {fmt.upper()} written to: {path}")
    print(f"CSV written to: {csv_path}")


if __name__ == "__main__":
    main()
