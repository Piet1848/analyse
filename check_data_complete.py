#!/usr/bin/env python3
"""
Check raw run folders for complete observable output.

The script walks a data root recursively, treating every directory with an
input.yaml as one run, matching the raw-run discovery used by search_data.py.
For each enabled observable it checks that the configured output file exists
and that the last data row starts with nSweep from MetropolisParams.
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Optional

from load_input_yaml import load_gradient_flow_params, load_params


@dataclass(frozen=True)
class FileCheck:
    name: str
    path: str
    expected_step: int
    last_step: Optional[int]
    status: str
    message: str


@dataclass(frozen=True)
class RunCheck:
    run_dir: str
    parent_dir: str
    run_name: str
    group_key: str
    params_label: str
    n_sweep: Optional[int]
    status: str
    checks: list[FileCheck]
    error: Optional[str] = None


def find_run_dirs(root: Path) -> list[Path]:
    return sorted(path.parent for path in root.rglob("input.yaml") if path.is_file())


def iter_expected_files(run_dir: Path, gauge: Any, gradient: Any) -> Iterable[tuple[str, Path]]:
    if gauge.measure_plaquette:
        yield "plaquette", run_dir / gauge.plaquette_filename
    if gauge.measure_wilson_loop_temporal:
        yield "wilson_temporal", run_dir / gauge.W_temp_filename
    if gauge.measure_wilson_loop_mu_nu:
        yield "wilson_mu_nu", run_dir / gauge.W_mu_nu_filename
    if gauge.measure_retrace_U:
        yield "retrace_U", run_dir / gauge.RetraceU_filename

    if gradient.enabled and gradient.measure_energy_clover:
        yield "gradient_flow", run_dir / gradient.obs_filename
    if gradient.enabled and gradient.measure_wilson_loop_temporal:
        yield "gradient_wilson_temporal", run_dir / gradient.W_temp_filename
    if gradient.enabled and gradient.measure_wilson_loop_mu_nu:
        yield "gradient_wilson_mu_nu", run_dir / gradient.W_mu_nu_filename
    if gradient.enabled and gradient.extract_t0:
        yield "gradient_t0", run_dir / gradient.t0_filename


def run_group_key(metro: Any, gauge: Any, gradient: Any) -> str:
    metro_dict = asdict(metro)
    metro_dict.pop("seed", None)
    payload = {
        "metro": metro_dict,
        "gauge": asdict(gauge),
        "gradient_flow": asdict(gradient),
    }
    return json.dumps(payload, sort_keys=True, separators=(",", ":"))


def params_label(metro: Any) -> str:
    lattice = "x".join(str(getattr(metro, f"L{i}")) for i in range(4))
    return (
        f"beta={metro.beta:g} L={lattice} "
        f"eps1={metro.epsilon1:g} eps2={metro.epsilon2:g} nSweep={metro.nSweep}"
    )


def read_last_step(path: Path) -> tuple[Optional[int], Optional[str]]:
    try:
        with path.open("rb") as handle:
            handle.seek(0, 2)
            position = handle.tell()
            buffer = b""
            while position > 0:
                chunk_size = min(8192, position)
                position -= chunk_size
                handle.seek(position)
                buffer = handle.read(chunk_size) + buffer
                for raw_line in reversed(buffer.splitlines()):
                    stripped = raw_line.strip()
                    if not stripped or stripped.startswith(b"#"):
                        continue
                    first_column = stripped.split(maxsplit=1)[0]
                    return int(float(first_column)), None
    except FileNotFoundError:
        return None, "missing file"
    except OSError as exc:
        return None, f"cannot read file: {exc}"
    except ValueError:
        return None, "last data row does not start with a numeric step"

    return None, "no data rows"


def check_file(name: str, path: Path, expected_step: int) -> FileCheck:
    last_step, error = read_last_step(path)
    if error:
        return FileCheck(
            name=name,
            path=str(path),
            expected_step=expected_step,
            last_step=last_step,
            status="missing" if error == "missing file" else "bad",
            message=error,
        )

    if last_step == expected_step:
        return FileCheck(
            name=name,
            path=str(path),
            expected_step=expected_step,
            last_step=last_step,
            status="ok",
            message="complete",
        )

    status = "short" if last_step is not None and last_step < expected_step else "over"
    return FileCheck(
        name=name,
        path=str(path),
        expected_step=expected_step,
        last_step=last_step,
        status=status,
        message=f"last step is {last_step}, expected {expected_step}",
    )


def check_run(run_dir: Path) -> RunCheck:
    try:
        metro, gauge = load_params(str(run_dir / "input.yaml"))
        gradient = load_gradient_flow_params(str(run_dir / "input.yaml"))
    except Exception as exc:
        return RunCheck(
            run_dir=str(run_dir),
            parent_dir=str(run_dir.parent),
            run_name=run_dir.name,
            group_key=f"error:{run_dir}",
            params_label="unreadable input.yaml",
            n_sweep=None,
            status="bad",
            checks=[],
            error=f"cannot load input.yaml: {exc}",
        )

    checks = [
        check_file(name, path, metro.nSweep)
        for name, path in iter_expected_files(run_dir, gauge, gradient)
    ]
    status = "ok" if checks and all(check.status == "ok" for check in checks) else "bad"
    if not checks:
        status = "empty"
    return RunCheck(
        run_dir=str(run_dir),
        parent_dir=str(run_dir.parent),
        run_name=run_dir.name,
        group_key=run_group_key(metro, gauge, gradient),
        params_label=params_label(metro),
        n_sweep=metro.nSweep,
        status=status,
        checks=checks,
    )


def natural_run_key(name: str) -> tuple[int, int | str]:
    try:
        return (0, int(name))
    except ValueError:
        return (1, name)


def compact_names(names: list[str]) -> str:
    numeric: list[int] = []
    text: list[str] = []
    for name in names:
        try:
            numeric.append(int(name))
        except ValueError:
            text.append(name)

    parts: list[str] = []
    sorted_numbers = sorted(set(numeric))
    if sorted_numbers:
        start = prev = sorted_numbers[0]
        for value in sorted_numbers[1:]:
            if value == prev + 1:
                prev = value
                continue
            parts.append(str(start) if start == prev else f"{start}-{prev}")
            start = prev = value
        parts.append(str(start) if start == prev else f"{start}-{prev}")

    parts.extend(sorted(text))
    return ",".join(parts)


def value_range(values: list[Optional[int]]) -> str:
    if any(value is None for value in values):
        return "None"
    ints = sorted({value for value in values if value is not None})
    if not ints:
        return "None"
    if len(ints) == 1:
        return str(ints[0])
    return f"{ints[0]}-{ints[-1]}"


def summarize_statuses(statuses: list[str]) -> str:
    counts: dict[str, int] = {}
    for status in statuses:
        counts[status] = counts.get(status, 0) + 1
    if len(counts) == 1:
        return next(iter(counts)).upper()
    return "MIXED " + ", ".join(f"{name}={count}" for name, count in sorted(counts.items()))


def print_grouped_report(results: list[RunCheck], only_incomplete: bool, show_details: bool) -> None:
    groups: dict[tuple[str, str, str], list[RunCheck]] = {}
    for run in results:
        if only_incomplete and run.status == "ok":
            continue
        groups.setdefault((run.parent_dir, run.group_key, run.status), []).append(run)

    if not groups:
        print("All checked runs are complete.")
        return

    for (parent_dir, _, status), runs in sorted(groups.items(), key=lambda item: item[0]):
        runs = sorted(runs, key=lambda run: natural_run_key(run.run_name))
        run_names = compact_names([run.run_name for run in runs])
        label = runs[0].params_label
        all_checks = [check for run in runs for check in run.checks]
        last = value_range([check.last_step for check in all_checks])
        expected = value_range([check.expected_step for check in all_checks])
        check_status = summarize_statuses([check.status for check in all_checks]) if all_checks else "NO FILES"
        print(
            f"{status.upper():7} {parent_dir}/{{{run_names}}}  runs={len(runs)}  "
            f"{label}  files={check_status} last={last} expected={expected}"
        )

        errors = [run for run in runs if run.error]
        for run in errors:
            print(f"  BAD     {run.run_name}: {run.error}")

        if not show_details:
            continue

        check_names = sorted({check.name for run in runs for check in run.checks})
        for check_name in check_names:
            checks = [check for run in runs for check in run.checks if check.name == check_name]
            statuses = [check.status for check in checks]
            expected = value_range([check.expected_step for check in checks])
            last = value_range([check.last_step for check in checks])
            print(
                f"  {summarize_statuses(statuses):<18} {check_name:<26} "
                f"last={last} expected={expected} files={len(checks)}/{len(runs)}"
            )


def print_text_report(results: list[RunCheck], only_incomplete: bool) -> None:
    shown = 0
    for run in results:
        if only_incomplete and run.status == "ok":
            continue
        shown += 1
        if run.error:
            print(f"{run.status.upper():7} {run.run_dir}: {run.error}")
            continue

        print(f"{run.status.upper():7} {run.run_dir}  nSweep={run.n_sweep}")
        for check in run.checks:
            if only_incomplete and check.status == "ok":
                continue
            print(
                f"  {check.status.upper():7} {check.name:<26} "
                f"last={check.last_step} expected={check.expected_step}  "
                f"{check.path}  {check.message}"
            )

    if shown == 0:
        print("All checked runs are complete.")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Check raw run folders for output files whose last step reaches nSweep."
    )
    parser.add_argument(
        "root",
        nargs="?",
        default="../data",
        help="Root data directory to scan recursively (default: ../data).",
    )
    parser.add_argument(
        "--only-incomplete",
        action="store_true",
        help="Show only runs/files that are missing, unreadable, or do not reach nSweep.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print machine-readable JSON instead of a text report.",
    )
    parser.add_argument(
        "--no-group",
        action="store_true",
        help="Print one entry per run instead of grouping equivalent runs.",
    )
    parser.add_argument(
        "--details",
        action="store_true",
        help="In grouped text output, also print the per-observable breakdown.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit with status 1 when any run is incomplete.",
    )
    args = parser.parse_args()

    root = Path(args.root).expanduser()
    if not root.exists():
        sys.exit(f"Root does not exist: {root}")

    results = [check_run(run_dir) for run_dir in find_run_dirs(root)]
    if args.json:
        print(json.dumps([asdict(result) for result in results], indent=2))
    elif args.no_group:
        print_text_report(results, args.only_incomplete)
    else:
        print_grouped_report(results, args.only_incomplete, args.details)

    if args.strict and any(result.status != "ok" for result in results):
        sys.exit(1)


if __name__ == "__main__":
    main()
