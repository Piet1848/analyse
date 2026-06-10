#!/usr/bin/env python3
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from analyze_precomputed_wilson_potential import DEFAULT_ANALYSIS_DIR, filename_token
from calculator import cornell_potential_ansatz, fit_r0_from_potential_data
from finalized_analysis_helpers import save_cornell_plot, save_json, save_r0_stability_plot


DEFAULT_FLOW_TIME = 0.5
DEFAULT_POTENTIAL_DIR = DEFAULT_ANALYSIS_DIR / f"potential_t_over_a2_{filename_token(DEFAULT_FLOW_TIME)}"


def parse_positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("value must be >= 1")
    return parsed


def parse_r_min_values(value: str) -> list[int]:
    values: set[int] = set()
    for chunk in value.split(","):
        item = chunk.strip()
        if not item:
            continue
        if ":" in item:
            parts = item.split(":")
            if len(parts) not in {2, 3}:
                raise argparse.ArgumentTypeError("ranges must be START:STOP or START:STOP:STEP")
            start = int(parts[0])
            stop = int(parts[1])
            step = int(parts[2]) if len(parts) == 3 else 1
            if step <= 0:
                raise argparse.ArgumentTypeError("range STEP must be positive")
            values.update(range(start, stop + 1, step))
        else:
            values.add(int(item))

    valid = sorted(value for value in values if value >= 1)
    if not valid:
        raise argparse.ArgumentTypeError("at least one r_min >= 1 is required")
    return valid


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        import json

        return json.load(handle)


def load_potential_points(potential_dir: Path) -> list[dict[str, Any]]:
    summary_path = potential_dir / "V_R_summary.json"
    if not summary_path.exists():
        raise FileNotFoundError(f"missing finalized potential summary: {summary_path}")

    payload = load_json(summary_path)
    results = payload.get("results", {})
    if not isinstance(results, dict) or not results:
        raise RuntimeError(f"no finalized V(R) points found in {summary_path}")

    points: list[dict[str, Any]] = []
    for r_key, row in results.items():
        try:
            r_value = int(row.get("R", r_key))
            value = float(row["value"])
        except (KeyError, TypeError, ValueError) as exc:
            raise RuntimeError(f"invalid V(R) row for R={r_key!r} in {summary_path}") from exc

        bootstrap_path = row.get("bootstrap_path")
        if bootstrap_path is not None:
            bootstrap_path = Path(str(bootstrap_path))
            if not bootstrap_path.is_absolute():
                bootstrap_path = potential_dir / bootstrap_path

        err = row.get("err")
        points.append(
            {
                "R": r_value,
                "V": value,
                "err": None if err is None else float(err),
                "flow_time": row.get("flow_time"),
                "fit_C": row.get("fit_C"),
                "fit_T": row.get("fit_T", []),
                "bootstrap_path": None if bootstrap_path is None else str(bootstrap_path),
            }
        )

    points = sorted(points, key=lambda item: int(item["R"]))
    if len(points) < 3:
        raise RuntimeError(f"need at least three finalized V(R) points, found {len(points)}")
    return points


def load_bootstrap_arrays(points: list[dict[str, Any]]) -> dict[int, np.ndarray]:
    arrays: dict[int, np.ndarray] = {}
    for point in points:
        path_text = point.get("bootstrap_path")
        if not path_text:
            continue
        path = Path(path_text)
        if not path.exists():
            continue
        arr = np.asarray(np.load(path), dtype=float)
        if arr.ndim != 1 or arr.size == 0:
            continue
        arrays[int(point["R"])] = arr
    return arrays


def selected_bootstrap_matrix(
    selected_points: list[dict[str, Any]],
    arrays_by_r: dict[int, np.ndarray],
    max_samples: int | None,
) -> np.ndarray | None:
    arrays = [arrays_by_r.get(int(point["R"])) for point in selected_points]
    if any(array is None for array in arrays):
        return None

    n_bootstrap = min(int(array.size) for array in arrays if array is not None)
    if max_samples is not None:
        n_bootstrap = min(n_bootstrap, int(max_samples))
    if n_bootstrap <= 0:
        return None

    return np.stack([array[:n_bootstrap] for array in arrays if array is not None], axis=0)


def finite_or_none(value: Any) -> float | None:
    if value is None:
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if np.isfinite(parsed) else None


def fit_for_r_min(
    points: list[dict[str, Any]],
    arrays_by_r: dict[int, np.ndarray],
    *,
    r_min: int,
    r_max: int | None,
    target_force: float,
    max_bootstrap_samples: int | None,
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]], np.ndarray]:
    selected = [
        point
        for point in points
        if int(point["R"]) >= int(r_min) and (r_max is None or int(point["R"]) <= int(r_max))
    ]
    if len(selected) < 3:
        raise ValueError(f"r_min={r_min} leaves only {len(selected)} V(R) point(s)")

    rs = np.asarray([float(point["R"]) for point in selected], dtype=float)
    vs = np.asarray([float(point["V"]) for point in selected], dtype=float)
    errs = np.asarray(
        [float(point["err"]) if point["err"] is not None else 1.0 for point in selected],
        dtype=float,
    )
    bootstrap_matrix = selected_bootstrap_matrix(selected, arrays_by_r, max_bootstrap_samples)
    fit = fit_r0_from_potential_data(
        rs,
        vs,
        errs=errs,
        bootstrap_matrix=bootstrap_matrix,
        target_force=target_force,
        n_bootstrap=0 if bootstrap_matrix is None else int(bootstrap_matrix.shape[1]),
    )

    params = fit["cornell_params"]
    r_values_all = [int(point["R"]) for point in points]
    fit_x = np.linspace(min(r_values_all), max(r_values_all), 300)
    curve_rows = [
        {
            "kind": "curve",
            "r_min": int(r_min),
            "R": float(r_value),
            "V_fit": float(cornell_potential_ansatz(r_value, params["A"], params["sigma"], params["B"])),
        }
        for r_value in fit_x
    ]
    point_rows = [
        {
            "kind": "point",
            "r_min": int(r_min),
            "R": int(point["R"]),
            "V": float(point["V"]),
            "err": point["err"],
            "used_in_fit": int(point["R"]) >= int(r_min) and (r_max is None or int(point["R"]) <= int(r_max)),
        }
        for point in points
    ]

    bootstrap_samples = np.asarray(fit.get("bootstrap_samples", []), dtype=float)
    fit_record = {
        "r_min": int(r_min),
        "r_max": None if r_max is None else int(r_max),
        "r0": float(fit["r0"]),
        "err": finite_or_none(fit.get("r0_err")),
        "r0_err": finite_or_none(fit.get("r0_err")),
        "target_force": float(target_force),
        "cornell_params": {
            "A": float(params["A"]),
            "sigma": float(params["sigma"]),
            "B": float(params["B"]),
        },
        "chi2": finite_or_none(fit.get("chi2")),
        "dof": int(fit["dof"]),
        "chi2_dof": finite_or_none(fit.get("chi2_dof")),
        "n_points": int(len(selected)),
        "fit_R": [int(value) for value in rs],
        "fit_V": [float(value) for value in vs],
        "fit_err": [float(value) for value in errs],
        "bootstrap_samples_used": int(bootstrap_matrix.shape[1]) if bootstrap_matrix is not None else 0,
    }
    return fit_record, point_rows, curve_rows, bootstrap_samples


def build_default_r_min_values(points: list[dict[str, Any]], r_max: int | None) -> list[int]:
    usable_rs = [int(point["R"]) for point in points if r_max is None or int(point["R"]) <= int(r_max)]
    return [r_value for idx, r_value in enumerate(usable_rs) if len(usable_rs[idx:]) >= 3]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fit the Cornell potential to finalized V(R) data from "
            "analyze_precomputed_wilson_potential.py, scanning r_min by default."
        )
    )
    parser.add_argument(
        "potential_dir",
        nargs="?",
        type=Path,
        default=DEFAULT_POTENTIAL_DIR,
        help=f"Directory containing V_R_summary.json. Default: {DEFAULT_POTENTIAL_DIR}",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Output directory. Default: <potential_dir>/cornell_fit.",
    )
    parser.add_argument(
        "--r-min",
        dest="r_min_values",
        type=parse_r_min_values,
        help="Comma-separated r_min values or inclusive ranges, e.g. 1,2,4:6. Default: all values leaving >=3 points.",
    )
    parser.add_argument("--r-max", type=parse_positive_int, help="Largest R included in every Cornell fit.")
    parser.add_argument(
        "--target-force",
        type=float,
        default=1.65,
        help="Sommer target in r0^2 F(r0). Default: 1.65.",
    )
    parser.add_argument(
        "--max-bootstrap-samples",
        type=parse_positive_int,
        help="Cap the number of V(R) bootstrap replicas used in each fit.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    potential_dir = args.potential_dir.resolve()
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else potential_dir / "cornell_fit"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    points = load_potential_points(potential_dir)
    arrays_by_r = load_bootstrap_arrays(points)

    r_min_values = (
        args.r_min_values
        if args.r_min_values is not None
        else build_default_r_min_values(points, args.r_max)
    )
    if not r_min_values:
        raise RuntimeError("no r_min values leave at least three V(R) points")

    r0_records: list[dict[str, Any]] = []
    curve_records: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    bootstrap_payload: dict[str, np.ndarray] = {}

    for r_min in r_min_values:
        try:
            fit_record, point_rows, curve_rows, bootstrap_samples = fit_for_r_min(
                points,
                arrays_by_r,
                r_min=int(r_min),
                r_max=args.r_max,
                target_force=args.target_force,
                max_bootstrap_samples=args.max_bootstrap_samples,
            )
        except (RuntimeError, ValueError) as exc:
            skipped.append({"r_min": int(r_min), "reason": str(exc)})
            continue

        r0_records.append(fit_record)
        curve_records.extend(point_rows)
        curve_records.extend(curve_rows)
        if bootstrap_samples.size:
            bootstrap_payload[f"r_min_{int(r_min)}"] = bootstrap_samples

    if not r0_records:
        details = "; ".join(f"r_min={row['r_min']}: {row['reason']}" for row in skipped)
        raise RuntimeError(f"no Cornell fits succeeded. {details}")

    r0_records = sorted(r0_records, key=lambda item: int(item["r_min"]))
    best_record = r0_records[0]
    result_payload = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source": "precomputed_wilson_potential",
        "potential_dir": str(potential_dir),
        "V_R_summary": str(potential_dir / "V_R_summary.json"),
        "target_force": float(args.target_force),
        "r_max": None if args.r_max is None else int(args.r_max),
        "n_potential_points": len(points),
        "potential_points": points,
        "fits": r0_records,
        "default_fit": best_record,
        "skipped": skipped,
        "outputs": {
            "results": str(output_dir / "cornell_fit_results.json"),
            "r0_stability_plot": str(output_dir / "r0_stability.html"),
            "cornell_plot": str(output_dir / "cornell_fit.html"),
        },
    }

    if bootstrap_payload:
        bootstrap_path = output_dir / "cornell_r0_bootstrap_samples.npz"
        np.savez(bootstrap_path, **bootstrap_payload)
        result_payload["outputs"]["r0_bootstrap_samples"] = str(bootstrap_path)

    save_json(output_dir / "cornell_fit_results.json", result_payload)
    save_json(output_dir / "cornell_curve_records.json", curve_records)
    save_r0_stability_plot(output_dir / "r0_stability.html", r0_records)
    save_cornell_plot(output_dir / "cornell_fit.html", curve_records, r0_records)

    print(f"Wrote Cornell fit analysis to {output_dir}")
    print("r_min  n  r0/a        err         sigma       B           chi2/dof")
    for row in r0_records:
        err_text = "nan" if row["err"] is None else f"{row['err']:.4g}"
        chi_text = "nan" if row["chi2_dof"] is None else f"{row['chi2_dof']:.4g}"
        print(
            f"{row['r_min']:>5d}  {row['n_points']:>1d}  "
            f"{row['r0']:<10.6g} {err_text:<11} "
            f"{row['cornell_params']['sigma']:<11.6g} "
            f"{row['cornell_params']['B']:<11.6g} {chi_text}"
        )
    if skipped:
        print(f"Skipped {len(skipped)} r_min value(s); see cornell_fit_results.json.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
