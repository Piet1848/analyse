from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from calculator import Calculator, cornell_potential_ansatz, fit_r0_from_potential_data


def json_default(value: Any):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def save_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, default=json_default)


def load_json(path: Path, default: Any = None) -> Any:
    if not path.exists():
        return default
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def summarize_bootstrap(values) -> tuple[float | None, float | None]:
    arr = np.asarray(values, dtype=float)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return None, None
    return float(np.mean(finite)), float(np.std(finite))


def enumerate_t_windows(unique_ts) -> list[tuple[int, int | None]]:
    t_values = sorted({int(t) for t in unique_ts})
    windows: list[tuple[int, int | None]] = []
    for idx, t_min in enumerate(t_values):
        remaining = t_values[idx:]
        if len(remaining) >= 3:
            windows.append((int(t_min), None))
        for t_max in remaining:
            usable = [t for t in t_values if int(t_min) <= t <= int(t_max)]
            if len(usable) >= 3:
                windows.append((int(t_min), int(t_max)))
    return sorted(set(windows), key=lambda item: (item[0], math.inf if item[1] is None else item[1]))


def window_label(t_min: int, t_max: int | None) -> str:
    return f"t_min={int(t_min)}, t_max={'None' if t_max is None else int(t_max)}"


def build_wrt_scan(calc: Calculator, unique_rs, unique_ts) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for r_val in unique_rs:
        for t_val in unique_ts:
            try:
                var = calc.get_variable("W_R_T", R=int(r_val), T=int(t_val))
            except Exception:
                continue
            value = var.get()
            if value is None or not np.isfinite(value):
                continue
            err = var.err()
            records.append(
                {
                    "R": int(r_val),
                    "T": int(t_val),
                    "value": float(value),
                    "err": float(err) if err is not None and np.isfinite(err) else None,
                }
            )
    return records


def build_effective_mass_scan(calc: Calculator, unique_rs, unique_ts) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    t_values = sorted({int(t) for t in unique_ts})
    if len(t_values) < 2:
        return records

    for r_val in unique_rs:
        for t_val in t_values[:-1]:
            try:
                var = calc.get_variable("effective_mass", R=int(r_val), T=int(t_val))
            except Exception:
                continue
            value = var.get()
            if value is None or not np.isfinite(value):
                continue
            err = var.err()
            records.append(
                {
                    "R": int(r_val),
                    "T": int(t_val),
                    "t_mid": float(t_val) + 0.5,
                    "value": float(value),
                    "err": float(err) if err is not None and np.isfinite(err) else None,
                }
            )
    return records


def build_v_r_scan(
    calc: Calculator,
    unique_rs,
    windows: list[tuple[int, int | None]],
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for t_min, t_max in windows:
        for r_val in unique_rs:
            try:
                var = calc.get_variable("V_R", R=int(r_val), t_min=int(t_min), t_max=t_max)
            except Exception:
                continue
            value = var.get()
            if value is None or not np.isfinite(value):
                continue
            err = var.err()
            fit_c = var.parameters.get("fit_C")
            records.append(
                {
                    "R": int(r_val),
                    "t_min": int(t_min),
                    "t_max": None if t_max is None else int(t_max),
                    "value": float(value),
                    "err": float(err) if err is not None and np.isfinite(err) else None,
                    "fit_C": float(fit_c) if fit_c is not None and np.isfinite(fit_c) else None,
                }
            )
    return records


def build_locked_r0_scan(
    v_results_by_r: dict[int, dict[str, Any]],
    target_force: float,
    n_bootstrap: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    r_values = sorted(int(r) for r in v_results_by_r)
    r0_records: list[dict[str, Any]] = []
    curve_records: list[dict[str, Any]] = []

    for r_min in r_values:
        fit_rs = [r for r in r_values if r >= int(r_min)]
        if len(fit_rs) < 3:
            continue

        rs = np.asarray(fit_rs, dtype=float)
        vs = np.asarray([v_results_by_r[r]["value"] for r in fit_rs], dtype=float)
        errs = np.asarray(
            [v_results_by_r[r].get("err") if v_results_by_r[r].get("err") is not None else 1.0 for r in fit_rs],
            dtype=float,
        )
        boot_matrix = np.stack(
            [np.asarray(v_results_by_r[r]["bootstrap_samples"], dtype=float) for r in fit_rs],
            axis=0,
        )

        try:
            fit_result = fit_r0_from_potential_data(
                rs,
                vs,
                errs=errs,
                bootstrap_matrix=boot_matrix,
                target_force=float(target_force),
                n_bootstrap=n_bootstrap,
            )
        except ValueError:
            continue

        r0_val = fit_result["r0"]
        if not np.isfinite(r0_val):
            continue

        params = fit_result["cornell_params"]
        r0_records.append(
            {
                "r_min": int(r_min),
                "r0": float(r0_val),
                "err": float(fit_result["r0_err"]) if fit_result["r0_err"] is not None else None,
                "A": float(params["A"]),
                "sigma": float(params["sigma"]),
                "B": float(params["B"]),
                "target_force": float(target_force),
            }
        )

        for r_val in fit_rs:
            curve_records.append(
                {
                    "kind": "point",
                    "r_min": int(r_min),
                    "R": float(r_val),
                    "V_fit": float(
                        cornell_potential_ansatz(
                            float(r_val),
                            float(params["A"]),
                            float(params["sigma"]),
                            float(params["B"]),
                        )
                    ),
                    "V": float(v_results_by_r[r_val]["value"]),
                    "err": (
                        float(v_results_by_r[r_val]["err"])
                        if v_results_by_r[r_val].get("err") is not None
                        else None
                    ),
                }
            )

        r_grid = np.linspace(float(min(fit_rs)), float(max(fit_rs)), 300)
        for r_val in r_grid:
            curve_records.append(
                {
                    "kind": "curve",
                    "r_min": int(r_min),
                    "R": float(r_val),
                    "V_fit": float(
                        cornell_potential_ansatz(
                            float(r_val),
                            float(params["A"]),
                            float(params["sigma"]),
                            float(params["B"]),
                        )
                    ),
                    "V": None,
                    "err": None,
                }
            )

    return r0_records, curve_records


def _get_plotly():
    try:
        import plotly.graph_objects as go
    except ImportError as exc:
        raise RuntimeError(
            "Plotly is required for finalize_analysis.py. Install it with `pip install plotly` "
            "or add it to your environment from requirements.txt."
        ) from exc
    return go


def _write_figure_html(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(path, include_plotlyjs="cdn", full_html=True)


def save_effective_mass_plot(
    path: Path,
    effective_mass_records: list[dict[str, Any]],
    v_scan_records: list[dict[str, Any]],
    r_value: int,
) -> None:
    go = _get_plotly()
    eff_rows = sorted(
        [row for row in effective_mass_records if int(row["R"]) == int(r_value)],
        key=lambda row: row["T"],
    )
    if not eff_rows:
        raise ValueError(f"No effective mass data available for R={r_value}")

    candidate_rows = sorted(
        [row for row in v_scan_records if int(row["R"]) == int(r_value)],
        key=lambda row: (row["t_min"], math.inf if row["t_max"] is None else row["t_max"]),
    )
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=[row["t_mid"] for row in eff_rows],
            y=[row["value"] for row in eff_rows],
            mode="lines+markers",
            name=f"m_eff (R={int(r_value)})",
            error_y={
                "type": "data",
                "array": [0.0 if row["err"] is None else row["err"] for row in eff_rows],
                "visible": any(row["err"] is not None for row in eff_rows),
            },
        )
    )

    window_to_trace_indices: dict[str, list[int]] = {}
    x_vals = [row["t_mid"] for row in eff_rows]
    x_min = min(x_vals) - 0.25
    x_max = max(x_vals) + 0.25

    if not candidate_rows:
        fig.update_layout(
            title=f"Effective mass for R={int(r_value)}",
            xaxis_title="T + 1/2",
            yaxis_title="m_eff",
            template="plotly_white",
        )
        _write_figure_html(fig, path)
        return

    for idx, row in enumerate(candidate_rows):
        label = window_label(int(row["t_min"]), row["t_max"])
        trace_indices = window_to_trace_indices.setdefault(label, [])
        v_val = float(row["value"])
        v_err = row["err"]
        visible = idx == 0

        trace_indices.append(len(fig.data))
        fig.add_trace(
            go.Scatter(
                x=[x_min, x_max],
                y=[v_val, v_val],
                mode="lines",
                name=f"V(R) {label}",
                visible=visible,
                line={"dash": "dash"},
            )
        )

        if v_err is not None and np.isfinite(v_err):
            trace_indices.append(len(fig.data))
            fig.add_trace(
                go.Scatter(
                    x=[x_min, x_max, x_max, x_min],
                    y=[v_val - v_err, v_val - v_err, v_val + v_err, v_val + v_err],
                    fill="toself",
                    mode="lines",
                    line={"width": 0},
                    fillcolor="rgba(31, 119, 180, 0.18)",
                    name=f"{label} error band",
                    visible=visible,
                    showlegend=False,
                )
            )

    base_visible = [True] + [False] * (len(fig.data) - 1)
    buttons = []
    for label, indices in window_to_trace_indices.items():
        visible = base_visible.copy()
        for trace_idx in indices:
            visible[trace_idx] = True
        buttons.append(
            {
                "label": label,
                "method": "update",
                "args": [
                    {"visible": visible},
                    {"title": f"Effective mass for R={int(r_value)} with {label}"},
                ],
            }
        )

    first_label = next(iter(window_to_trace_indices))
    fig.update_layout(
        title=f"Effective mass for R={int(r_value)} with {first_label}",
        xaxis_title="T + 1/2",
        yaxis_title="m_eff",
        template="plotly_white",
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
        margin={"l": 60, "r": 240, "t": 70, "b": 60},
    )
    _write_figure_html(fig, path)


def save_r0_stability_plot(path: Path, r0_records: list[dict[str, Any]]) -> None:
    go = _get_plotly()
    rows = sorted(r0_records, key=lambda row: row["r_min"])
    if not rows:
        raise ValueError("No r0 scan data available")

    fig = go.Figure(
        data=[
            go.Scatter(
                x=[row["r_min"] for row in rows],
                y=[row["r0"] for row in rows],
                mode="lines+markers",
                error_y={
                    "type": "data",
                    "array": [0.0 if row["err"] is None else row["err"] for row in rows],
                    "visible": any(row["err"] is not None for row in rows),
                },
                name="r0/a",
            )
        ]
    )
    fig.update_layout(
        title="r0 stability by r_min",
        xaxis_title="r_min",
        yaxis_title="r0/a",
        template="plotly_white",
    )
    _write_figure_html(fig, path)


def save_cornell_plot(
    path: Path,
    curve_records: list[dict[str, Any]],
    r0_records: list[dict[str, Any]],
) -> None:
    go = _get_plotly()
    if not curve_records or not r0_records:
        raise ValueError("No Cornell fit data available")

    r_min_values = [int(row["r_min"]) for row in sorted(r0_records, key=lambda row: row["r_min"])]
    fig = go.Figure()
    trace_indices: dict[int, list[int]] = {}

    for idx, r_min in enumerate(r_min_values):
        visible = idx == 0
        trace_indices[r_min] = []
        point_rows = sorted(
            [row for row in curve_records if row["kind"] == "point" and int(row["r_min"]) == r_min],
            key=lambda row: row["R"],
        )
        curve_rows = sorted(
            [row for row in curve_records if row["kind"] == "curve" and int(row["r_min"]) == r_min],
            key=lambda row: row["R"],
        )
        if curve_rows:
            trace_indices[r_min].append(len(fig.data))
            fig.add_trace(
                go.Scatter(
                    x=[row["R"] for row in curve_rows],
                    y=[row["V_fit"] for row in curve_rows],
                    mode="lines",
                    name=f"Cornell fit (r_min={r_min})",
                    visible=visible,
                )
            )
        if point_rows:
            trace_indices[r_min].append(len(fig.data))
            fig.add_trace(
                go.Scatter(
                    x=[row["R"] for row in point_rows],
                    y=[row["V"] for row in point_rows],
                    mode="markers",
                    name=f"Selected V(R) (r_min={r_min})",
                    visible=visible,
                    error_y={
                        "type": "data",
                        "array": [0.0 if row["err"] is None else row["err"] for row in point_rows],
                        "visible": any(row["err"] is not None for row in point_rows),
                    },
                )
            )

    buttons = []
    for row in sorted(r0_records, key=lambda item: item["r_min"]):
        r_min = int(row["r_min"])
        visible = [False] * len(fig.data)
        for trace_idx in trace_indices.get(r_min, []):
            visible[trace_idx] = True
        buttons.append(
            {
                "label": f"r_min={r_min}",
                "method": "update",
                "args": [
                    {"visible": visible},
                    {"title": f"Cornell fit for r_min={r_min} (r0/a={row['r0']:.4f})"},
                ],
            }
        )

    first = sorted(r0_records, key=lambda item: item["r_min"])[0]
    fig.update_layout(
        title=f"Cornell fit for r_min={int(first['r_min'])} (r0/a={first['r0']:.4f})",
        xaxis_title="R",
        yaxis_title="V(R)",
        template="plotly_white",
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
        margin={"l": 60, "r": 240, "t": 70, "b": 60},
    )
    _write_figure_html(fig, path)
