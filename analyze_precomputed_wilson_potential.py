#!/usr/bin/env python3
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
import warnings
import webbrowser

import numpy as np
from scipy.optimize import curve_fit

from finalized_analysis_helpers import (
    _get_plotly,
    _write_figure_html,
    save_effective_mass_plot,
    save_json,
)


DEFAULT_ANALYSIS_DIR = (
    Path(__file__).resolve().parent.parent
    / "data"
    / "gradient_flow_wtemp_analysis"
    / "beta_2p4__L0_24__epsilon1_0__dt_0p01__b6c467be"
)


def filename_token(value: float | int | str) -> str:
    return str(value).replace("-", "m").replace("+", "").replace(".", "p")


def exponential_ansatz(t, c_0, v):
    exponent = np.clip(-float(v) * np.asarray(t, dtype=float), -745.0, 700.0)
    with np.errstate(over="ignore", invalid="ignore"):
        return float(c_0) * np.exp(exponent)


def parse_t_min(value: str) -> int:
    t_min = int(value)
    if t_min < 1:
        raise argparse.ArgumentTypeError("t_min must be >= 1")
    return t_min


def read_precomputed_wilson_stats(analysis_dir: Path, flow_time: float) -> list[dict[str, Any]]:
    stats_path = analysis_dir / "wilson_loop_stats.dat"
    if not stats_path.exists():
        raise FileNotFoundError(f"missing precomputed stats file: {stats_path}")

    rows: list[dict[str, Any]] = []
    with stats_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) < 12:
                continue

            row_flow = float(parts[0])
            if not np.isclose(row_flow, float(flow_time)):
                continue

            bootstrap_path = Path(parts[11])
            if not bootstrap_path.is_absolute():
                bootstrap_path = analysis_dir / bootstrap_path

            rows.append(
                {
                    "t_over_a2": row_flow,
                    "R": int(parts[1]),
                    "T": int(parts[2]),
                    "n_measurements": int(parts[3]),
                    "first_conf_id": int(float(parts[4])),
                    "last_conf_id": int(float(parts[5])),
                    "block_size": int(parts[6]),
                    "n_blocks": int(parts[7]),
                    "tail": int(parts[8]),
                    "value": float(parts[9]),
                    "err": float(parts[10]),
                    "bootstrap_path": str(bootstrap_path),
                    "flow_time": float(flow_time),
                }
            )

    if not rows:
        raise RuntimeError(f"no Wilson-loop rows found for t_over_a2={flow_time:g} in {stats_path}")
    return sorted(rows, key=lambda row: (int(row["R"]), int(row["T"])))


def build_effective_mass_records(wrt_records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_pair = {(int(row["R"]), int(row["T"])): row for row in wrt_records}
    records: list[dict[str, Any]] = []

    for r_value in sorted({int(row["R"]) for row in wrt_records}):
        t_values = sorted(int(row["T"]) for row in wrt_records if int(row["R"]) == r_value)
        for t_value in t_values:
            current = by_pair.get((r_value, t_value))
            following = by_pair.get((r_value, t_value + 1))
            if current is None or following is None:
                continue
            if current["value"] <= 0 or following["value"] <= 0:
                continue

            current_boot = np.asarray(np.load(current["bootstrap_path"]), dtype=float)
            following_boot = np.asarray(np.load(following["bootstrap_path"]), dtype=float)
            n_boot = min(current_boot.size, following_boot.size)
            current_boot = current_boot[:n_boot]
            following_boot = following_boot[:n_boot]

            boot = np.full(n_boot, np.nan, dtype=float)
            valid = (
                np.isfinite(current_boot)
                & np.isfinite(following_boot)
                & (current_boot > 0)
                & (following_boot > 0)
            )
            boot[valid] = np.log(current_boot[valid] / following_boot[valid])
            finite_boot = boot[np.isfinite(boot)]

            records.append(
                {
                    "R": r_value,
                    "T": t_value,
                    "t_mid": float(t_value) + 0.5,
                    "value": float(np.log(current["value"] / following["value"])),
                    "err": float(np.std(finite_boot)) if finite_boot.size else None,
                    "flow_time": float(current["flow_time"]),
                    "bootstrap_finite_fraction": float(finite_boot.size / boot.size) if boot.size else None,
                }
            )

    return records


def _initial_exponential_guess(ts: np.ndarray, ws: np.ndarray) -> tuple[float, float]:
    try:
        p = np.polyfit(ts, np.log(ws), 1)
        return max(float(np.exp(p[1])), 1e-300), float(-p[0])
    except (RuntimeError, ValueError, TypeError, np.linalg.LinAlgError, FloatingPointError):
        c_guess = max(float(ws[0]), 1e-300)
        v_guess = -float(np.log(ws[-1] / ws[0]) / (ts[-1] - ts[0])) if ts[-1] != ts[0] else 0.1
        return c_guess, v_guess


def _fit_exponential(
    ts: np.ndarray,
    ws: np.ndarray,
    sigma: np.ndarray | None = None,
    p0: tuple[float, float] | None = None,
) -> tuple[float, float]:
    if len(ws) == 0 or np.any(~np.isfinite(ws)) or np.any(ws <= 0):
        return np.nan, np.nan

    try:
        c_guess, v_guess = p0 if p0 is not None else _initial_exponential_guess(ts, ws)
        fit_sigma = sigma if sigma is not None and np.all(np.isfinite(sigma)) and np.all(sigma > 0) else None
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            popt, _ = curve_fit(
                exponential_ansatz,
                ts,
                ws,
                p0=[max(float(c_guess), 1e-300), float(v_guess)],
                sigma=fit_sigma,
                absolute_sigma=(fit_sigma is not None),
                bounds=([0.0, -np.inf], [np.inf, np.inf]),
                maxfev=5000,
            )
        c_fit, v_fit = popt
        return float(v_fit), float(c_fit)
    except (RuntimeError, ValueError, TypeError, np.linalg.LinAlgError):
        return np.nan, np.nan


def fit_potential_window(
    rows_for_r: list[dict[str, Any]],
    t_min: int,
    t_max: int | None,
    *,
    include_bootstrap: bool = True,
) -> dict[str, Any] | None:
    selected = [
        row
        for row in sorted(rows_for_r, key=lambda item: int(item["T"]))
        if int(row["T"]) >= int(t_min) and (t_max is None or int(row["T"]) <= int(t_max))
    ]
    if len(selected) < 3:
        return None

    ts_fit = np.asarray([int(row["T"]) for row in selected], dtype=float)
    ws_main = np.asarray([float(row["value"]) for row in selected], dtype=float)
    ws_err = np.asarray([float(row["err"]) for row in selected], dtype=float)

    bad_indices = np.where((~np.isfinite(ws_main)) | (ws_main <= 0))[0]
    if len(bad_indices) > 0:
        first_bad = int(bad_indices[0])
        ts_fit = ts_fit[:first_bad]
        ws_main = ws_main[:first_bad]
        ws_err = ws_err[:first_bad]
        selected = selected[:first_bad]

    if len(ts_fit) < 3:
        return None

    if np.any((~np.isfinite(ws_err)) | (ws_err <= 0)):
        positive = ws_err[np.isfinite(ws_err) & (ws_err > 0)]
        replacement = float(np.mean(positive)) if positive.size else 1.0
        ws_err = ws_err.copy()
        ws_err[(~np.isfinite(ws_err)) | (ws_err <= 0)] = replacement

    v_main, c_main = _fit_exponential(ts_fit, ws_main, sigma=ws_err)
    if not np.isfinite(v_main):
        return None

    result: dict[str, Any] = {
        "R": int(selected[0]["R"]),
        "t_min": int(t_min),
        "t_max": None if t_max is None else int(t_max),
        "value": float(v_main),
        "err": None,
        "fit_C": float(c_main),
        "flow_time": float(selected[0]["flow_time"]),
        "fit_T": [int(value) for value in ts_fit],
    }
    if not include_bootstrap:
        return result

    boot_arrays = [np.asarray(np.load(row["bootstrap_path"]), dtype=float) for row in selected]
    n_bootstrap = min(array.size for array in boot_arrays)
    if n_bootstrap <= 0:
        return None

    all_boots = np.stack([array[:n_bootstrap] for array in boot_arrays], axis=0)
    bootstrap_vs = np.empty(n_bootstrap, dtype=float)
    bootstrap_p0 = (c_main, v_main)

    for idx in range(n_bootstrap):
        ws_sample = np.array(all_boots[:, idx], copy=True)
        invalid = (~np.isfinite(ws_sample)) | (ws_sample <= 0)
        if np.any(invalid):
            ws_sample[invalid] = ws_main[invalid]
        if np.any(ws_sample <= 0):
            bootstrap_vs[idx] = v_main
            continue

        v_sample, _ = _fit_exponential(ts_fit, ws_sample, sigma=ws_err, p0=bootstrap_p0)
        bootstrap_vs[idx] = v_sample if np.isfinite(v_sample) else v_main

    finite_boot = bootstrap_vs[np.isfinite(bootstrap_vs)]
    result["err"] = float(np.std(finite_boot)) if finite_boot.size else None
    result["bootstrap_samples"] = bootstrap_vs
    return result


def build_v_r_scan(wrt_records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for r_value in sorted({int(row["R"]) for row in wrt_records}):
        rows_for_r = [row for row in wrt_records if int(row["R"]) == r_value]
        t_values = sorted({int(row["T"]) for row in rows_for_r})
        for idx, t_min in enumerate(t_values):
            if len(t_values[idx:]) < 3:
                continue
            result = fit_potential_window(rows_for_r, t_min, None, include_bootstrap=True)
            if result is None:
                continue
            record = dict(result)
            record.pop("bootstrap_samples", None)
            records.append(record)
    return records


def save_v_r_window_scan_plot(path: Path, scan_records: list[dict[str, Any]]) -> None:
    go = _get_plotly()
    r_values = sorted({int(row["R"]) for row in scan_records})
    if not r_values:
        raise ValueError("No V(R) candidate data available")

    fig = go.Figure()
    trace_indices_by_r: dict[int, list[int]] = {}
    for idx, r_value in enumerate(r_values):
        visible = idx == 0
        rows_for_r = sorted(
            [row for row in scan_records if int(row["R"]) == r_value],
            key=lambda row: int(row["t_min"]),
        )
        trace_indices_by_r[r_value] = []
        if rows_for_r:
            trace_indices_by_r[r_value].append(len(fig.data))
            fig.add_trace(
                go.Scatter(
                    x=[int(row["t_min"]) for row in rows_for_r],
                    y=[float(row["value"]) for row in rows_for_r],
                    mode="lines+markers",
                    name=f"R={r_value}",
                    visible=visible,
                    error_y={
                        "type": "data",
                        "array": [0.0 if row["err"] is None else float(row["err"]) for row in rows_for_r],
                        "visible": any(row["err"] is not None for row in rows_for_r),
                    },
                    hovertemplate=(
                        "R=%{customdata[0]}<br>"
                        "t_min=%{x}<br>"
                        "t_max=None<br>"
                        "V=%{y:.8g}<br>"
                        "err=%{customdata[1]:.4g}<extra></extra>"
                    ),
                    customdata=[
                        [r_value, 0.0 if row["err"] is None else float(row["err"])]
                        for row in rows_for_r
                    ],
                )
            )

    buttons = []
    for r_value in r_values:
        visible = [False] * len(fig.data)
        for trace_idx in trace_indices_by_r[r_value]:
            visible[trace_idx] = True
        buttons.append(
            {
                "label": f"R={r_value}",
                "method": "update",
                "args": [
                    {"visible": visible},
                    {"title": f"V(R) window scan for R={r_value}"},
                ],
            }
        )

    first_r = r_values[0]
    fig.update_layout(
        title=f"V(R) window scan for R={first_r}",
        xaxis_title="t_min",
        yaxis_title="V(R)",
        template="plotly_white",
        margin={"l": 60, "r": 180, "t": 80, "b": 60},
        updatemenus=[
            {
                "buttons": buttons,
                "direction": "down",
                "showactive": True,
                "x": 1.02,
                "xanchor": "left",
                "y": 1.0,
                "yanchor": "top",
            }
        ],
    )
    fig.update_xaxes(dtick=1)
    _write_figure_html(fig, path)


def open_plot(path: Path) -> None:
    uri = path.resolve().as_uri()
    try:
        webbrowser.open(uri)
    except Exception as exc:
        print(f"Could not open plot automatically: {exc}")
    print(f"Plot: {path}")


def _prompt_t_min(r_value: int, candidates: list[dict[str, Any]], previous: int | None) -> int | None:
    if not candidates:
        print(f"\nR = {r_value}: no finite candidate windows")
        return None

    print(f"\nR = {r_value}")
    while True:
        prompt = "Enter t_min"
        if previous is not None:
            prompt += f", or press Enter to reuse t_min={previous}"
        prompt += "; type skip to skip this R: "
        response = input(prompt).strip()
        if response.lower() == "skip":
            return None
        if not response and previous is not None:
            return previous
        try:
            choice = parse_t_min(response)
        except (ValueError, argparse.ArgumentTypeError):
            print("Please enter an integer t_min.")
            continue
        valid = {int(row["t_min"]) for row in candidates}
        if choice in valid:
            return choice
        print(f"t_min={choice} is not available for R={r_value}.")


def finalize_potentials(
    wrt_records: list[dict[str, Any]],
    scan_records: list[dict[str, Any]],
    output_dir: Path,
    *,
    t_min: int | None = None,
    interactive: bool = False,
) -> dict[str, Any]:
    selected_choices: dict[str, Any] = {"choices": {}}
    summary: dict[str, Any] = {"results": {}}
    bootstrap_dir = output_dir / "bootstrap"
    previous: int | None = None

    for r_value in sorted({int(row["R"]) for row in wrt_records}):
        rows_for_r = [row for row in wrt_records if int(row["R"]) == r_value]
        candidates = [row for row in scan_records if int(row["R"]) == r_value]

        selected_t_min = t_min
        if selected_t_min is None and interactive:
            selected_t_min = _prompt_t_min(r_value, candidates, previous)
        if selected_t_min is None:
            continue

        candidate_t_mins = {int(row["t_min"]) for row in candidates}
        if selected_t_min not in candidate_t_mins:
            print(f"Skipping R={r_value}: t_min={selected_t_min} has no finite candidate fit.")
            continue

        result = fit_potential_window(rows_for_r, selected_t_min, None)
        if result is None:
            print(f"Skipping R={r_value}: selected window did not produce a finite fit.")
            continue

        previous = selected_t_min
        bootstrap_dir.mkdir(parents=True, exist_ok=True)
        bootstrap_path = bootstrap_dir / f"V_R_R_{r_value}.npy"
        np.save(bootstrap_path, np.asarray(result["bootstrap_samples"], dtype=float))

        selected_choices["choices"][str(r_value)] = {
            "R": r_value,
            "t_min": int(selected_t_min),
            "t_max": None,
            "flow_time": float(result["flow_time"]),
            "saved_at": datetime.now(timezone.utc).isoformat(),
        }
        summary["results"][str(r_value)] = {
            "value": float(result["value"]),
            "err": result["err"],
            "fit_C": float(result["fit_C"]),
            "flow_time": float(result["flow_time"]),
            "fit_T": result["fit_T"],
            "bootstrap_path": str(bootstrap_path),
        }

    save_json(output_dir / "choices.json", selected_choices)
    save_json(output_dir / "V_R_summary.json", summary)
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract V(R) from precomputed gradient_flow_wtemp_analysis output. "
            "Uses effective masses only as diagnostics and fits W(R,T)=C exp[-V T] directly."
        )
    )
    parser.add_argument(
        "analysis_dir",
        nargs="?",
        type=Path,
        default=DEFAULT_ANALYSIS_DIR,
        help=f"Directory containing wilson_loop_stats.dat. Default: {DEFAULT_ANALYSIS_DIR}",
    )
    parser.add_argument("--flow-time", type=float, default=0.5, help="Wilson-loop t_over_a2 value to analyze.")
    parser.add_argument("--output-dir", type=Path, help="Output directory. Default: <analysis_dir>/potential_t_over_a2_<flow>.")
    parser.add_argument("--t-min", type=parse_t_min, help="Finalize all R with this t_min, fitting all T >= t_min.")
    parser.add_argument("--interactive", action="store_true", help="Open plots and prompt for t_min for each R.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    analysis_dir = args.analysis_dir.resolve()
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else analysis_dir / f"potential_t_over_a2_{filename_token(args.flow_time)}"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    wrt_records = read_precomputed_wilson_stats(analysis_dir, args.flow_time)
    effective_mass_records = build_effective_mass_records(wrt_records)
    scan_records = build_v_r_scan(wrt_records)

    save_json(output_dir / "wrt_scan.json", wrt_records)
    save_json(output_dir / "effective_mass_scan.json", effective_mass_records)
    save_json(output_dir / "v_r_scan_candidates.json", scan_records)
    effective_mass_plot = output_dir / "effective_mass_all_R.html"
    window_scan_plot = output_dir / "v_r_window_scan.html"
    save_effective_mass_plot(effective_mass_plot, wrt_records, effective_mass_records)
    save_v_r_window_scan_plot(window_scan_plot, scan_records)

    if args.interactive:
        print("Use these plots to choose the fit windows:")
        open_plot(effective_mass_plot)
        open_plot(window_scan_plot)

    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source": "precomputed_gradient_flow_wtemp_analysis",
        "analysis_dir": str(analysis_dir),
        "flow_time": float(args.flow_time),
        "n_wilson_loop_records": len(wrt_records),
        "n_effective_mass_records": len(effective_mass_records),
        "n_v_r_candidates": len(scan_records),
        "outputs": {
            "wrt_scan": str(output_dir / "wrt_scan.json"),
            "effective_mass_scan": str(output_dir / "effective_mass_scan.json"),
            "v_r_scan_candidates": str(output_dir / "v_r_scan_candidates.json"),
            "effective_mass_plot": str(effective_mass_plot),
            "window_scan_plot": str(window_scan_plot),
        },
    }

    if args.t_min is not None or args.interactive:
        summary = finalize_potentials(
            wrt_records,
            scan_records,
            output_dir,
            t_min=args.t_min,
            interactive=args.interactive,
        )
        manifest["outputs"]["choices"] = str(output_dir / "choices.json")
        manifest["outputs"]["V_R_summary"] = str(output_dir / "V_R_summary.json")
        manifest["n_final_v_r"] = len(summary.get("results", {}))

    save_json(output_dir / "manifest.json", manifest)

    print(f"Wrote precomputed Wilson-potential analysis to {output_dir}")
    if not args.interactive:
        print(f"Effective-mass diagnostic plot: {effective_mass_plot}")
        print(f"Window-scan plot: {window_scan_plot}")
    if "n_final_v_r" in manifest:
        print(f"Finalized {manifest['n_final_v_r']} V(R) point(s).")
    else:
        print("No V(R) points finalized. Use --t-min or --interactive to write V_R_summary.json.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
