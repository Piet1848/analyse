from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

import data_organizer as do
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


def build_wrt_scan(calc: Calculator, unique_rs, unique_ts, flow_time: float | None = None) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for r_val in unique_rs:
        for t_val in unique_ts:
            try:
                var = calc.get_variable("W_R_T", R=int(r_val), T=int(t_val), flow_time=flow_time)
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
                    "flow_time": None if flow_time is None else float(flow_time),
                }
            )
    return records


def build_bootstrap_block_size_scan(
    file_data,
    block_sizes: list[int],
    *,
    n_bootstrap: int,
    pair_keys: list[tuple[int, int]] | None = None,
    flow_time: float | None = None,
) -> tuple[list[dict[str, Any]], list[tuple[int, int]]]:
    probe_calc = Calculator(file_data, n_bootstrap=1, step_size=1)
    available_pairs = probe_calc.get_available_pairs(flow_time=flow_time)
    selected_pairs = (
        _select_preview_pair_keys([(int(r_val), int(t_val)) for r_val, t_val in available_pairs])
        if pair_keys is None
        else [(int(r_val), int(t_val)) for r_val, t_val in pair_keys]
    )
    if not selected_pairs:
        raise ValueError("No W(R,T) pairs available for bootstrap block-size scan.")

    records: list[dict[str, Any]] = []
    for block_size in sorted({int(value) for value in block_sizes if int(value) >= 1}):
        calc = Calculator(file_data, n_bootstrap=int(n_bootstrap), step_size=int(block_size))
        try:
            calc.prime_w_rt_cache(pairs=selected_pairs, flow_time=flow_time)
        except Exception:
            pass

        for r_val, t_val in selected_pairs:
            row: dict[str, Any] = {
                "block_size": int(block_size),
                "R": int(r_val),
                "T": int(t_val),
                "value": None,
                "err": None,
                "status": "ok",
                "flow_time": None if flow_time is None else float(flow_time),
            }
            try:
                var = calc.get_variable("W_R_T", R=int(r_val), T=int(t_val), flow_time=flow_time)
                value = var.get()
                err = var.err()
                row["value"] = float(value) if value is not None and np.isfinite(value) else None
                row["err"] = float(err) if err is not None and np.isfinite(err) else None
            except Exception as exc:
                row["status"] = str(exc)
            records.append(row)

    return records, selected_pairs


def recommend_bootstrap_block_size(
    scan_records: list[dict[str, Any]],
    *,
    fallback: int,
    relative_tolerance: float = 0.10,
) -> int:
    by_pair: dict[tuple[int, int], dict[int, float]] = {}
    for row in scan_records:
        err = row.get("err")
        if err is None or not np.isfinite(err):
            continue
        by_pair.setdefault((int(row["R"]), int(row["T"])), {})[int(row["block_size"])] = float(err)

    block_sizes = sorted({block for pair_errors in by_pair.values() for block in pair_errors})
    if len(block_sizes) < 2:
        return max(1, int(fallback))

    for current, next_block in zip(block_sizes, block_sizes[1:], strict=False):
        relative_changes: list[float] = []
        for pair_errors in by_pair.values():
            current_err = pair_errors.get(current)
            next_err = pair_errors.get(next_block)
            if current_err is None or next_err is None or current_err <= 0:
                continue
            relative_changes.append(abs(next_err - current_err) / current_err)
        if relative_changes and float(np.median(relative_changes)) <= float(relative_tolerance):
            return int(current)

    return max(1, int(fallback))


def build_thermalization_preview(
    run_dirs: list[str],
    *,
    include_combined: bool = True,
) -> list[dict[str, Any]]:
    run_series: list[dict[str, Any]] = []

    for run_dir in run_dirs:
        print(f"Processing run directory for thermalization preview: {run_dir}")
        w_temp_path = Path(run_dir) / "W_temp.out"
        if not w_temp_path.exists():
            continue

        file_data = do.FileData(str(w_temp_path))
        file_data.read_file()
        file_data.align_lengths()

        if not file_data.observables:
            continue
        if min((len(obs.values) for obs in file_data.observables), default=0) <= 0:
            continue

        try:
            obs_w = np.asarray(file_data.get("W_temp").values, dtype=float)
            obs_l = np.asarray(file_data.get("L").values)
            obs_t = np.asarray(file_data.get("T").values)
        except ValueError:
            continue

        step_obs = next((o for o in file_data.observables if o.name in {"# step", "step"}), None)
        if step_obs is None:
            continue

        steps = np.asarray(step_obs.values, dtype=float)
        if steps.size == 0:
            continue

        rows_per_cfg = do._infer_rows_per_configuration(steps, obs_l, obs_t)
        if rows_per_cfg <= 0:
            continue
        if len(obs_w) % rows_per_cfg != 0 or len(steps) % rows_per_cfg != 0:
            continue

        pair_order = [(int(l_val), int(t_val)) for l_val, t_val in zip(obs_l[:rows_per_cfg], obs_t[:rows_per_cfg])]
        if len(set(pair_order)) != len(pair_order):
            continue

        step_series = steps[::rows_per_cfg]
        w_matrix = obs_w.reshape(-1, rows_per_cfg)
        run_label = Path(run_dir).name

        for idx, (r_val, t_val) in enumerate(pair_order):
            values = np.asarray(w_matrix[:, idx], dtype=float)
            if values.size == 0:
                continue
            running_mean = np.cumsum(values) / np.arange(1, values.size + 1, dtype=float)
            run_series.append(
                {
                    "kind": "run",
                    "run_label": run_label,
                    "R": int(r_val),
                    "T": int(t_val),
                    "steps": [float(step) for step in step_series],
                    "values": [float(value) for value in values],
                    "running_mean": [float(value) for value in running_mean],
                }
            )
    print(f"Collected thermalization preview data for {len(run_series)} series from {len(run_dirs)} run directories.")
    if not include_combined:
        return run_series

    combined_series: list[dict[str, Any]] = []
    pair_keys = sorted({(int(row["R"]), int(row["T"])) for row in run_series}, key=lambda item: (item[0], item[1]))
    for r_val, t_val in pair_keys:
        print(f"Constructing combined series for R={r_val}, T={t_val} from {sum(1 for row in run_series if int(row['R']) == r_val and int(row['T']) == t_val)} runs.")
        grouped: dict[float, list[float]] = {}
        for row in run_series:
            if int(row["R"]) != r_val or int(row["T"]) != t_val:
                continue
            for step, value in zip(row["steps"], row["values"]):
                grouped.setdefault(float(step), []).append(float(value))

        if not grouped:
            continue

        sorted_steps = sorted(grouped)
        means = np.asarray([np.mean(grouped[step]) for step in sorted_steps], dtype=float)
        errs = np.asarray([np.std(grouped[step]) for step in sorted_steps], dtype=float)
        running_mean = np.cumsum(means) / np.arange(1, len(means) + 1, dtype=float)
        combined_series.append(
            {
                "kind": "combined",
                "run_label": "combined",
                "R": int(r_val),
                "T": int(t_val),
                "steps": [float(step) for step in sorted_steps],
                "values": [float(value) for value in means],
                "err": [float(value) for value in errs],
                "running_mean": [float(value) for value in running_mean],
            }
        )
    print(f"Constructed combined thermalization series for {len(combined_series)} R-T pairs.")
    return combined_series + run_series


def build_effective_mass_scan(calc: Calculator, unique_rs, unique_ts, flow_time: float | None = None) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    t_values = sorted({int(t) for t in unique_ts})
    if len(t_values) < 2:
        return records

    for r_val in unique_rs:
        for t_val in t_values[:-1]:
            try:
                var = calc.get_variable("effective_mass", R=int(r_val), T=int(t_val), flow_time=flow_time)
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
                    "flow_time": None if flow_time is None else float(flow_time),
                }
            )
    return records


def build_v_r_scan(
    calc: Calculator,
    unique_rs,
    windows: list[tuple[int, int | None]],
    flow_time: float | None = None,
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for t_min, t_max in windows:
        for r_val in unique_rs:
            try:
                var = calc.get_variable(
                    "V_R",
                    R=int(r_val),
                    t_min=int(t_min),
                    t_max=t_max,
                    flow_time=flow_time,
                )
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
                    "flow_time": None if flow_time is None else float(flow_time),
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
                "chi2": float(fit_result["chi2"]) if fit_result["chi2"] is not None else None,
                "dof": int(fit_result["dof"]),
                "chi2_dof": float(fit_result["chi2_dof"]) if fit_result["chi2_dof"] is not None else None,
            }
        )

        for r_val in r_values:
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
                    "used_in_fit": bool(r_val in fit_rs),
                }
            )

        r_grid_min = min(1.0, float(min(fit_rs)))
        r_grid = np.linspace(r_grid_min, float(max(fit_rs)), 300)
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


def _select_preview_pair_keys(
    pair_keys: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    ordered = sorted(pair_keys, key=lambda item: (item[0], item[1]))
    if not ordered:
        return []

    available = set(ordered)
    max_r = max(r_val for r_val, _ in ordered)
    max_t_for_max_r = max(t_val for r_val, t_val in ordered if r_val == max_r)
    preferred = [(1, 1), (1, 8), (5, 5), (5, 8), (8, 10), (10, 10), (max_r, max_t_for_max_r)]

    selected: list[tuple[int, int]] = []
    for pair in preferred:
        if pair in available and pair not in selected:
            selected.append(pair)

    if selected:
        return selected

    if len(ordered) <= 5:
        return ordered

    selected_indices: list[int] = []
    for idx in range(5):
        position = round(idx * (len(ordered) - 1) / 4)
        if position not in selected_indices:
            selected_indices.append(position)

    for position in range(len(ordered)):
        if len(selected_indices) >= 5:
            break
        if position not in selected_indices:
            selected_indices.append(position)

    return [ordered[position] for position in sorted(selected_indices[:5])]


def _trim_preview_series(
    row: dict[str, Any],
    *,
    fraction: float = 0.10,
    min_points: int = 3,
) -> dict[str, Any]:
    steps = [float(step) for step in row["steps"]]
    if not steps:
        return row

    keep_count = min(len(steps), max(min_points, int(math.ceil(len(steps) * float(fraction)))))
    trimmed = dict(row)
    trimmed["steps"] = steps[:keep_count]
    trimmed["values"] = [float(value) for value in row["values"][:keep_count]]
    trimmed["running_mean"] = [float(value) for value in row["running_mean"][:keep_count]]
    if "err" in row:
        trimmed["err"] = [float(value) for value in row.get("err", [])[:keep_count]]
    return trimmed


def save_thermalization_plot(
    path: Path,
    thermalization_records: list[dict[str, Any]],
    suggested_cut: int | None,
) -> None:
    go = _get_plotly()
    from plotly.subplots import make_subplots

    all_pair_keys = sorted(
        {(int(row["R"]), int(row["T"])) for row in thermalization_records},
        key=lambda item: (item[0], item[1]),
    )
    if not all_pair_keys:
        raise ValueError("No thermalization preview data available")
    pair_keys = _select_preview_pair_keys(all_pair_keys)
    selected_pairs = set(pair_keys)
    preview_records = [
        _trim_preview_series(row, fraction=0.10)
        for row in thermalization_records
        if (int(row["R"]), int(row["T"])) in selected_pairs
    ]

    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.09,
        subplot_titles=("Raw W(R,T) vs Monte Carlo step", "Running mean of W(R,T)"),
    )

    trace_indices: dict[tuple[int, int], list[int]] = {}

    for idx, pair_key in enumerate(pair_keys):
        r_val, t_val = pair_key
        visible = idx == 0
        trace_indices[pair_key] = []

        pair_rows = sorted(
            [row for row in preview_records if int(row["R"]) == r_val and int(row["T"]) == t_val],
            key=lambda row: (0 if row["kind"] == "combined" else 1, row["run_label"]),
        )

        for row in pair_rows:
            kind = str(row["kind"])
            run_label = str(row["run_label"])
            steps = [float(step) for step in row["steps"]]
            values = [float(value) for value in row["values"]]
            running_mean = [float(value) for value in row["running_mean"]]
            legend_name = "Combined mean" if kind == "combined" else f"Run {run_label}"
            line_width = 3 if kind == "combined" else 1.4
            opacity = 1.0 if kind == "combined" else 0.45
            legend_group = f"{kind}:{run_label}"

            if kind == "combined":
                err = [float(value) for value in row.get("err", [])]
                if err:
                    lower = [val - delta for val, delta in zip(values, err)]
                    upper = [val + delta for val, delta in zip(values, err)]
                    trace_indices[pair_key].append(len(fig.data))
                    fig.add_trace(
                        go.Scatter(
                            x=steps + steps[::-1],
                            y=lower + upper[::-1],
                            fill="toself",
                            mode="lines",
                            line={"width": 0},
                            fillcolor="rgba(31, 119, 180, 0.12)",
                            name="Combined spread",
                            hoverinfo="skip",
                            showlegend=False,
                            visible=visible,
                            legendgroup=legend_group,
                        ),
                        row=1,
                        col=1,
                    )

            trace_indices[pair_key].append(len(fig.data))
            fig.add_trace(
                go.Scatter(
                    x=steps,
                    y=values,
                    mode="lines+markers",
                    name=legend_name,
                    visible=visible,
                    opacity=opacity,
                    line={"width": line_width},
                    legendgroup=legend_group,
                    showlegend=True,
                ),
                row=1,
                col=1,
            )

            trace_indices[pair_key].append(len(fig.data))
            fig.add_trace(
                go.Scatter(
                    x=steps,
                    y=running_mean,
                    mode="lines",
                    name=f"{legend_name} running mean",
                    visible=visible,
                    opacity=opacity,
                    line={"width": line_width},
                    legendgroup=legend_group,
                    showlegend=False,
                ),
                row=2,
                col=1,
            )

    buttons = []
    for r_val, t_val in pair_keys:
        visible = [False] * len(fig.data)
        for trace_idx in trace_indices[(r_val, t_val)]:
            visible[trace_idx] = True
        buttons.append(
            {
                "label": f"R={r_val}, T={t_val}",
                "method": "update",
                "args": [
                    {"visible": visible},
                    {"title": f"Thermalization preview for W(R={r_val}, T={t_val})"},
                ],
            }
        )

    first_r, first_t = pair_keys[0]
    layout_updates: dict[str, Any] = {
        "title": f"Thermalization preview for W(R={first_r}, T={first_t})",
        "template": "plotly_white",
        "xaxis2_title": "Monte Carlo step",
        "yaxis_title": "W(R,T)",
        "yaxis2_title": "Running mean",
        "margin": {"l": 65, "r": 250, "t": 80, "b": 65},
        "annotations": [
            {
                "x": 0.0,
                "y": 1.12,
                "xref": "paper",
                "yref": "paper",
                "text": (
                    f"Preview limited to {len(pair_keys)} representative (R,T) pairs "
                    f"and the first 10% of Monte Carlo steps."
                ),
                "showarrow": False,
                "xanchor": "left",
                "font": {"size": 12, "color": "#555"},
            }
        ],
        "updatemenus": [
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
    }
    if suggested_cut is not None:
        layout_updates["shapes"] = [
            {
                "type": "line",
                "xref": "x",
                "yref": "paper",
                "x0": float(suggested_cut),
                "x1": float(suggested_cut),
                "y0": 0.0,
                "y1": 1.0,
                "line": {"color": "firebrick", "dash": "dash"},
            }
        ]
        layout_updates["annotations"].append(
            {
                "x": float(suggested_cut),
                "y": 1.02,
                "xref": "x",
                "yref": "paper",
                "text": f"suggested cut = {int(suggested_cut)}",
                "showarrow": False,
                "font": {"color": "firebrick"},
            }
        )

    fig.update_layout(
        **layout_updates,
    )
    _write_figure_html(fig, path)


def save_bootstrap_block_size_plot(
    path: Path,
    scan_records: list[dict[str, Any]],
    *,
    recommended_block_size: int | None,
) -> None:
    go = _get_plotly()
    finite_rows = [
        row
        for row in scan_records
        if row.get("err") is not None and np.isfinite(row["err"])
    ]
    if not finite_rows:
        raise ValueError("No finite bootstrap error data available for block-size plot.")

    pair_keys = sorted({(int(row["R"]), int(row["T"])) for row in finite_rows})
    fig = go.Figure()
    for r_val, t_val in pair_keys:
        rows = sorted(
            [
                row
                for row in finite_rows
                if int(row["R"]) == int(r_val) and int(row["T"]) == int(t_val)
            ],
            key=lambda row: int(row["block_size"]),
        )
        fig.add_trace(
            go.Scatter(
                x=[int(row["block_size"]) for row in rows],
                y=[float(row["err"]) for row in rows],
                mode="lines+markers",
                name=f"W(R={r_val}, T={t_val})",
            )
        )

    shapes: list[dict[str, Any]] = []
    annotations: list[dict[str, Any]] = [
        {
            "x": 0.0,
            "y": 1.12,
            "xref": "paper",
            "yref": "paper",
            "text": "Block size is measured in saved configurations, not raw Monte Carlo sweeps.",
            "showarrow": False,
            "xanchor": "left",
            "font": {"size": 12, "color": "#555"},
        }
    ]
    if recommended_block_size is not None:
        shapes.append(
            {
                "type": "line",
                "xref": "x",
                "yref": "paper",
                "x0": float(recommended_block_size),
                "x1": float(recommended_block_size),
                "y0": 0.0,
                "y1": 1.0,
                "line": {"color": "firebrick", "dash": "dash"},
            }
        )
        annotations.append(
            {
                "x": float(recommended_block_size),
                "y": 1.02,
                "xref": "x",
                "yref": "paper",
                "text": f"recommended block size = {int(recommended_block_size)}",
                "showarrow": False,
                "font": {"color": "firebrick"},
            }
        )

    fig.update_layout(
        title="Bootstrap error vs block size",
        xaxis_title="Bootstrap block size (saved configurations)",
        yaxis_title="Bootstrap error of W(R,T)",
        template="plotly_white",
        shapes=shapes,
        annotations=annotations,
        margin={"l": 70, "r": 180, "t": 90, "b": 70},
    )
    fig.update_xaxes(type="linear")
    _write_figure_html(fig, path)


def save_gradient_flow_plot(path: Path, flow_summary: dict[str, Any]) -> None:
    go = _get_plotly()
    times = [float(value) for value in flow_summary.get("available_flow_times", [])]
    t2e = flow_summary.get("t2E_clover", {})
    t2e_err = flow_summary.get("t2E_clover_err", {})
    values = [t2e.get(f"{time:.12g}") for time in times]
    errors = [t2e_err.get(f"{time:.12g}") for time in times]

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=times,
            y=values,
            error_y={"type": "data", "array": errors, "visible": True},
            mode="lines+markers",
            name="t2E_clover",
        )
    )
    target = flow_summary.get("t0_target")
    if target is not None:
        fig.add_hline(y=float(target), line_dash="dash", line_color="firebrick")
    if target is None or not np.isclose(float(target), 0.1):
        fig.add_hline(y=0.1, line_dash="dash", line_color="darkorange")
    t0 = flow_summary.get("t0")
    if t0 is not None:
        fig.add_vline(x=float(t0), line_dash="dot", line_color="seagreen")
    t2e_0p1 = flow_summary.get("t_over_a2_at_t2E_clover_0p1")
    if t2e_0p1 is not None:
        fig.add_vline(x=float(t2e_0p1), line_dash="dot", line_color="darkorange")
    t2e_0p1_fit = flow_summary.get("t_over_a2_at_t2E_clover_0p1_weighted_fit")
    if t2e_0p1_fit is not None:
        fig.add_vline(x=float(t2e_0p1_fit), line_dash="dashdot", line_color="purple")
    fig.update_layout(
        title="Gradient-flow energy observable",
        template="plotly_white",
        xaxis_title="t/a^2",
        yaxis_title="t^2 E_clover",
    )
    _write_figure_html(fig, path)


def save_creutz_plot(path: Path, creutz_summary: dict[str, Any]) -> None:
    go = _get_plotly()
    rows = _creutz_rows(creutz_summary)

    fig = go.Figure()
    if rows:
        rows = sorted(rows)
        fig.add_trace(
            go.Scatter(
                x=[row[0] for row in rows],
                y=[row[2] for row in rows],
                error_y={
                    "type": "data",
                    "array": [row[3] for row in rows],
                    "visible": True,
                },
                mode="markers",
                marker={"color": [row[1] for row in rows], "colorscale": "Viridis", "showscale": True},
                name="chi",
                text=[f"T={row[1]:g}" for row in rows],
            )
        )
    else:
        fig.add_annotation(
            text=creutz_summary.get("status", "No Creutz-ratio data available"),
            x=0.5,
            y=0.5,
            xref="paper",
            yref="paper",
            showarrow=False,
        )
    fig.update_layout(
        title="Creutz ratios",
        template="plotly_white",
        xaxis_title="R + 1/2",
        yaxis_title="chi",
    )
    _write_figure_html(fig, path)


def _creutz_rows(creutz_summary: dict[str, Any]) -> list[tuple[float, float, float, float | None]]:
    chi = creutz_summary.get("chi", {})
    chi_err = creutz_summary.get("chi_err", {})
    rows = []
    for key, value in chi.items():
        try:
            r_text, t_text = key.split(",", 1)
            err = chi_err.get(key)
            rows.append(
                (
                    float(r_text),
                    float(t_text),
                    float(value),
                    float(err) if err is not None else None,
                )
            )
        except (TypeError, ValueError):
            continue
    return rows


def save_creutz_diagonal_plot(path: Path, creutz_summary: dict[str, Any]) -> None:
    go = _get_plotly()
    rows = [
        row
        for row in _creutz_rows(creutz_summary)
        if np.isclose(row[0], row[1])
    ]

    fig = go.Figure()
    if rows:
        rows = sorted(rows)
        fig.add_trace(
            go.Scatter(
                x=[row[0] for row in rows],
                y=[row[2] for row in rows],
                error_y={
                    "type": "data",
                    "array": [row[3] for row in rows],
                    "visible": True,
                },
                mode="lines+markers",
                name="chi(R,R)",
            )
        )
    else:
        fig.add_annotation(
            text="No diagonal Creutz-ratio data available",
            x=0.5,
            y=0.5,
            xref="paper",
            yref="paper",
            showarrow=False,
        )
    fig.update_layout(
        title="Creutz ratios on R = T",
        template="plotly_white",
        xaxis_title="R = T midpoint",
        yaxis_title="chi(R,T)",
    )
    _write_figure_html(fig, path)


def save_effective_mass_plot(
    path: Path,
    wrt_records: list[dict[str, Any]],
    effective_mass_records: list[dict[str, Any]],
) -> None:
    go = _get_plotly()
    from plotly.subplots import make_subplots

    r_values = sorted({int(row["R"]) for row in effective_mass_records})
    if not r_values:
        raise ValueError("No effective mass data available")

    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=False,
        vertical_spacing=0.09,
        subplot_titles=("Effective mass plateau", "log W(R,T) used in the fit"),
    )

    trace_indices_by_r: dict[int, list[int]] = {}
    for idx, r_value in enumerate(r_values):
        visible = idx == 0
        trace_indices_by_r[int(r_value)] = []

        eff_rows = sorted(
            [row for row in effective_mass_records if int(row["R"]) == int(r_value)],
            key=lambda row: row["T"],
        )
        if eff_rows:
            trace_indices_by_r[int(r_value)].append(len(fig.data))
            fig.add_trace(
                go.Scatter(
                    x=[row["t_mid"] for row in eff_rows],
                    y=[row["value"] for row in eff_rows],
                    mode="lines+markers",
                    name=f"m_eff (R={int(r_value)})",
                    visible=visible,
                    error_y={
                        "type": "data",
                        "array": [0.0 if row["err"] is None else row["err"] for row in eff_rows],
                        "visible": any(row["err"] is not None for row in eff_rows),
                    },
                ),
                row=1,
                col=1,
            )

        wrt_rows = sorted(
            [row for row in wrt_records if int(row["R"]) == int(r_value)],
            key=lambda row: row["T"],
        )
        positive_wrt_rows = [
            row for row in wrt_rows
            if np.isfinite(row["value"]) and float(row["value"]) > 0
        ]
        if positive_wrt_rows:
            log_w_values = [float(np.log(row["value"])) for row in positive_wrt_rows]
            log_w_errors = [
                (float(row["err"]) / float(row["value"]))
                if row.get("err") is not None and np.isfinite(row["err"]) and float(row["err"]) >= 0
                else 0.0
                for row in positive_wrt_rows
            ]
            trace_indices_by_r[int(r_value)].append(len(fig.data))
            fig.add_trace(
                go.Scatter(
                    x=[row["T"] for row in positive_wrt_rows],
                    y=log_w_values,
                    mode="lines+markers",
                    name=f"log W(R={int(r_value)},T)",
                    visible=visible,
                    error_y={
                        "type": "data",
                        "array": log_w_errors,
                        "visible": any(err > 0 for err in log_w_errors),
                    },
                ),
                row=2,
                col=1,
            )

    buttons = []
    for r_value in r_values:
        visible = [False] * len(fig.data)
        for trace_idx in trace_indices_by_r[int(r_value)]:
            visible[trace_idx] = True
        buttons.append(
            {
                "label": f"R={int(r_value)}",
                "method": "update",
                "args": [
                    {"visible": visible},
                    {"title": f"Effective mass for R={int(r_value)}"},
                ],
            }
        )

    first_r = int(r_values[0])
    fig.update_layout(
        title=f"Effective mass for R={first_r}",
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
    fig.update_xaxes(title_text="T + 1/2", row=1, col=1)
    fig.update_yaxes(title_text="m_eff", row=1, col=1)
    fig.update_xaxes(title_text="T", row=2, col=1)
    fig.update_yaxes(title_text="log W(R,T)", row=2, col=1)
    _write_figure_html(fig, path)


def save_r0_stability_plot(path: Path, r0_records: list[dict[str, Any]]) -> None:
    go = _get_plotly()
    rows = sorted(r0_records, key=lambda row: row["r_min"])
    if not rows:
        raise ValueError("No r0 scan data available")
    r_min_values = [int(row["r_min"]) for row in rows]

    fig = go.Figure(
        data=[
            go.Scatter(
                x=r_min_values,
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
    fig.update_xaxes(
        tickmode="array",
        tickvals=r_min_values,
        range=[min(r_min_values) - 0.25, max(r_min_values) + 0.25],
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

    def format_title(row: dict[str, Any]) -> str:
        title = f"Cornell fit for r_min={int(row['r_min'])} (r0/a={row['r0']:.4f})"
        chi2_dof = row.get("chi2_dof")
        if chi2_dof is not None and np.isfinite(chi2_dof):
            title += f", chi2/dof={chi2_dof:.3f}"
        elif row.get("chi2") is not None and np.isfinite(row["chi2"]):
            title += f", chi2={float(row['chi2']):.3f}"
        return title

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
            used_rows = [row for row in point_rows if bool(row.get("used_in_fit"))]
            excluded_rows = [row for row in point_rows if not bool(row.get("used_in_fit"))]

            if excluded_rows:
                trace_indices[r_min].append(len(fig.data))
                fig.add_trace(
                    go.Scatter(
                        x=[row["R"] for row in excluded_rows],
                        y=[row["V"] for row in excluded_rows],
                        mode="markers",
                        name=f"Excluded V(R) (r_min={r_min})",
                        visible=visible,
                        marker={"symbol": "circle-open", "size": 9},
                        error_y={
                            "type": "data",
                            "array": [0.0 if row["err"] is None else row["err"] for row in excluded_rows],
                            "visible": any(row["err"] is not None for row in excluded_rows),
                        },
                    )
                )
            if used_rows:
                trace_indices[r_min].append(len(fig.data))
                fig.add_trace(
                    go.Scatter(
                        x=[row["R"] for row in used_rows],
                        y=[row["V"] for row in used_rows],
                        mode="markers",
                        name=f"Fit points V(R) (r_min={r_min})",
                        visible=visible,
                        marker={"size": 9},
                        error_y={
                            "type": "data",
                            "array": [0.0 if row["err"] is None else row["err"] for row in used_rows],
                            "visible": any(row["err"] is not None for row in used_rows),
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
                    {"title": format_title(row)},
                ],
            }
        )

    first = sorted(r0_records, key=lambda item: item["r_min"])[0]
    fig.update_layout(
        title=format_title(first),
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
