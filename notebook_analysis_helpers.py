from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from typing import Any
from datetime import datetime, timezone

import numpy as np

import run_evaluation
from calculator import Calculator, get_key_name, get_key_params
from finalized_analysis_helpers import (
    build_effective_mass_scan,
    build_locked_r0_scan,
    build_v_r_scan,
    build_wrt_scan,
    enumerate_t_windows,
    load_json,
    save_json,
)


CACHE_SCHEMA_VERSION = 1


def _canonical_run_dirs(run_dirs: list[str]) -> list[str]:
    return [os.path.abspath(path) for path in run_dirs]


def _json_hash(payload: Any) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()[:16]


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def input_hash_for_runs(run_dirs: list[str]) -> str:
    payload: list[dict[str, Any]] = []
    for run_dir in _canonical_run_dirs(run_dirs):
        input_path = Path(run_dir) / "input.yaml"
        entry: dict[str, Any] = {"run_dir": run_dir, "input_yaml_exists": input_path.exists()}
        if input_path.exists():
            entry["input_yaml_sha256"] = hashlib.sha256(input_path.read_bytes()).hexdigest()
        payload.append(entry)
    return _json_hash(payload)


def resolve_thermalization_steps_by_run(
    run_dirs: list[str],
    *,
    thermalization_steps: int | None = None,
    thermalization_steps_by_run: dict[str, int] | None = None,
) -> dict[str, int]:
    cuts = {
        os.path.abspath(path): int(value)
        for path, value in (thermalization_steps_by_run or {}).items()
    }
    default_cut = (
        int(run_evaluation.THERMALIZATION_STEPS)
        if thermalization_steps is None
        else int(thermalization_steps)
    )
    resolved: dict[str, int] = {}
    for run_dir in _canonical_run_dirs(run_dirs):
        cut = int(cuts.get(run_dir, default_cut))
        if cut < 0:
            raise ValueError(f"Thermalization cut must be non-negative for run: {run_dir}")
        resolved[run_dir] = cut
    return resolved


def build_analysis_options(
    *,
    thermalization_steps_by_run: dict[str, int],
    n_bootstrap: int = run_evaluation.DEFAULT_N_BOOTSTRAP,
    block_size: int = run_evaluation.DEFAULT_BLOCK_SIZE,
    target_force: float = run_evaluation.DEFAULT_SOMMER_TARGET,
    r0_t_min: int = run_evaluation.DEFAULT_V_R_T_MIN,
    r0_t_max: int | None = run_evaluation.DEFAULT_V_R_T_MAX,
    load_workers: int = 1,
    calc_workers: int | None = None,
    selected_v_windows: dict[int | str, tuple[int, int | None]] | None = None,
) -> dict[str, Any]:
    normalized_windows = None
    if selected_v_windows:
        normalized_windows = {
            str(int(r_value)): {
                "t_min": int(window[0]),
                "t_max": None if window[1] is None else int(window[1]),
            }
            for r_value, window in sorted(selected_v_windows.items(), key=lambda item: int(item[0]))
        }

    return {
        "thermalization_steps_by_run": {
            os.path.abspath(path): int(cut)
            for path, cut in sorted(thermalization_steps_by_run.items())
        },
        "n_bootstrap": int(n_bootstrap),
        "block_size": int(block_size),
        "target_force": float(target_force),
        "r0_t_min": int(r0_t_min),
        "r0_t_max": None if r0_t_max is None else int(r0_t_max),
        "load_workers": int(load_workers),
        "calc_workers": None if calc_workers is None else int(calc_workers),
        "selected_v_windows": normalized_windows,
    }


def build_scan_cache_key(
    run_dirs: list[str],
    analysis_options: dict[str, Any],
    *,
    input_hash: str | None = None,
) -> str:
    return _json_hash(
        {
            "source_run_dirs": _canonical_run_dirs(run_dirs),
            "input_hash": input_hash or input_hash_for_runs(run_dirs),
            "analysis_options": analysis_options,
            "calculator_version": run_evaluation.CALC_VERSION,
            "cache_schema_version": CACHE_SCHEMA_VERSION,
        }
    )


def default_scan_cache_path(
    run_dirs: list[str],
    analysis_options: dict[str, Any],
    *,
    cache_root: Path | None = None,
    input_hash: str | None = None,
) -> Path:
    root = run_evaluation.CALC_RESULT_BASE if cache_root is None else Path(cache_root)
    key = build_scan_cache_key(run_dirs, analysis_options, input_hash=input_hash)
    return root / f"notebook_scan_{key}.json"


def expected_cache_metadata(
    run_dirs: list[str],
    analysis_options: dict[str, Any],
    *,
    input_hash: str | None = None,
) -> dict[str, Any]:
    return {
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "source_run_dirs": _canonical_run_dirs(run_dirs),
        "input_hash": input_hash or input_hash_for_runs(run_dirs),
        "analysis_options": analysis_options,
        "calculator_version": run_evaluation.CALC_VERSION,
    }


def explain_scan_cache_staleness(
    cache_payload: dict[str, Any] | None,
    expected_metadata: dict[str, Any],
) -> list[str]:
    if cache_payload is None:
        return ["cache file is missing"]
    reasons: list[str] = []
    for field, expected_value in expected_metadata.items():
        actual_value = cache_payload.get(field)
        if actual_value != expected_value:
            reasons.append(f"{field} differs")
    if not isinstance(cache_payload.get("wrt_scan"), list):
        reasons.append("wrt_scan is missing")
    if not isinstance(cache_payload.get("effective_mass_scan"), list):
        reasons.append("effective_mass_scan is missing")
    if not isinstance(cache_payload.get("v_r_scan"), list):
        reasons.append("v_r_scan is missing")
    return reasons


def _collect_bootstrap_quality(calc: Calculator) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for key, var in sorted(calc.variables.items(), key=lambda item: repr(item[0])):
        boot = var.bootstrap()
        if boot is None:
            continue
        boot_arr = np.asarray(boot, dtype=float)
        rows.append(
            {
                "variable": get_key_name(key),
                "params": get_key_params(key),
                "n_bootstrap": int(boot_arr.size),
                "finite_fraction": float(np.isfinite(boot_arr).sum() / boot_arr.size) if boot_arr.size else None,
                "invalid_count_before_repair": int(getattr(var, "bootstrap_invalid_count", 0)),
                "error": float(var.err()) if var.err() is not None and np.isfinite(var.err()) else None,
            }
        )
    return rows


def _selected_v_results(
    calc: Calculator,
    unique_rs: list[int],
    *,
    t_min: int,
    t_max: int | None,
    selected_v_windows: dict[int | str, tuple[int, int | None]] | None = None,
) -> dict[int, dict[str, Any]]:
    results: dict[int, dict[str, Any]] = {}
    for r_value in unique_rs:
        selected = selected_v_windows.get(r_value) if selected_v_windows else None
        if selected is None and selected_v_windows:
            selected = selected_v_windows.get(str(r_value))
        current_t_min, current_t_max = selected if selected is not None else (t_min, t_max)
        try:
            var = calc.get_variable(
                "V_R",
                R=int(r_value),
                t_min=int(current_t_min),
                t_max=current_t_max,
            )
        except Exception:
            continue
        value = var.get()
        boot = var.bootstrap()
        if value is None or not np.isfinite(value) or boot is None:
            continue
        results[int(r_value)] = {
            "R": int(r_value),
            "t_min": int(current_t_min),
            "t_max": None if current_t_max is None else int(current_t_max),
            "value": float(value),
            "err": float(var.err()) if var.err() is not None and np.isfinite(var.err()) else None,
            "bootstrap_samples": np.asarray(boot, dtype=float),
        }
    return results


def build_physics_audit(cache_payload: dict[str, Any]) -> dict[str, Any]:
    wrt_records = list(cache_payload.get("wrt_scan", []))
    v_r_records = list(cache_payload.get("v_r_scan", []))
    r0_records = list(cache_payload.get("r0_scan", []))
    bootstrap_quality = list(cache_payload.get("bootstrap_quality", []))
    analysis_options = dict(cache_payload.get("analysis_options", {}))
    target_force = float(analysis_options.get("target_force", run_evaluation.DEFAULT_SOMMER_TARGET))

    audit: dict[str, Any] = {
        "summary": {"warning_count": 0},
        "cornell_sommer": {"warnings": [], "checks": []},
        "fit_weights": {"warnings": [], "checks": []},
        "creutz_ratio": {"warnings": [], "checks": []},
        "bootstrap_quality": {"warnings": [], "checks": []},
        "derived_units": {"warnings": [], "checks": []},
    }

    for row in r0_records:
        sigma = float(row.get("sigma", np.nan))
        b_param = float(row.get("B", np.nan))
        r0_val = float(row.get("r0", np.nan))
        check = {
            "r_min": row.get("r_min"),
            "sigma_positive": bool(np.isfinite(sigma) and sigma > 0),
            "B_less_than_target": bool(np.isfinite(b_param) and b_param < target_force),
            "chi2_dof": row.get("chi2_dof"),
        }
        if np.isfinite(sigma) and np.isfinite(b_param) and np.isfinite(r0_val):
            check["r0_squared_force"] = float((r0_val ** 2) * (sigma + b_param / (r0_val ** 2)))
            check["target_force_delta"] = float(check["r0_squared_force"] - target_force)
        if not check["sigma_positive"]:
            audit["cornell_sommer"]["warnings"].append(f"r_min={row.get('r_min')} has non-positive sigma")
        if not check["B_less_than_target"]:
            audit["cornell_sommer"]["warnings"].append(f"r_min={row.get('r_min')} has B >= target_force")
        if row.get("chi2_dof") is None:
            audit["cornell_sommer"]["warnings"].append(f"r_min={row.get('r_min')} has no chi2/dof")
        audit["cornell_sommer"]["checks"].append(check)
    if not r0_records:
        audit["cornell_sommer"]["warnings"].append("No r0 scan records available for Cornell/Sommer audit")

    rows_by_r: dict[int, list[dict[str, Any]]] = {}
    for row in wrt_records:
        rows_by_r.setdefault(int(row["R"]), []).append(row)
    for r_value, rows in sorted(rows_by_r.items()):
        non_positive = sum(1 for row in rows if float(row.get("value", np.nan)) <= 0)
        zero_or_missing_err = sum(
            1
            for row in rows
            if row.get("err") is None or not np.isfinite(row["err"]) or float(row["err"]) <= 0
        )
        check = {
            "R": int(r_value),
            "points": len(rows),
            "non_positive_wilson_points": non_positive,
            "zero_or_missing_error_points": zero_or_missing_err,
            "log_weight_convention": "np.polyfit receives inverse log-space std-dev weights W/err(W)",
        }
        if non_positive:
            audit["fit_weights"]["warnings"].append(f"R={r_value} has non-positive W(R,T) points")
        if zero_or_missing_err:
            audit["fit_weights"]["warnings"].append(f"R={r_value} has zero or missing W(R,T) errors")
        audit["fit_weights"]["checks"].append(check)

    wrt_lookup = {
        (int(row["R"]), int(row["T"])): float(row["value"])
        for row in wrt_records
        if row.get("value") is not None and np.isfinite(row["value"])
    }
    creutz_values: list[float] = []
    for r_value, t_value in sorted(wrt_lookup):
        required = [
            (r_value, t_value),
            (r_value + 1, t_value),
            (r_value, t_value + 1),
            (r_value + 1, t_value + 1),
        ]
        if not all(pair in wrt_lookup for pair in required):
            continue
        w_rt = wrt_lookup[(r_value, t_value)]
        w_r1t = wrt_lookup[(r_value + 1, t_value)]
        w_rt1 = wrt_lookup[(r_value, t_value + 1)]
        w_r1t1 = wrt_lookup[(r_value + 1, t_value + 1)]
        numerator = w_rt1 * w_r1t
        denominator = w_rt * w_r1t1
        if numerator > 0 and denominator > 0:
            creutz_values.append(float(np.log(numerator / denominator)))
    if creutz_values:
        median_chi = float(np.median(creutz_values))
        audit["creutz_ratio"]["checks"].append(
            {
                "n_adjacent_rectangles": len(creutz_values),
                "median_current_chi": median_chi,
                "convention": "log(W(R,T+1)*W(R+1,T)/(W(R,T)*W(R+1,T+1)))",
            }
        )
        if median_chi <= 0:
            audit["creutz_ratio"]["warnings"].append("Median Creutz ratio is non-positive")
    else:
        audit["creutz_ratio"]["warnings"].append("No adjacent Wilson-loop rectangles available for Creutz audit")

    repaired_total = sum(int(row.get("invalid_count_before_repair") or 0) for row in bootstrap_quality)
    low_finite = [
        row
        for row in bootstrap_quality
        if row.get("finite_fraction") is not None and float(row["finite_fraction"]) < 0.95
    ]
    audit["bootstrap_quality"]["checks"].append(
        {
            "variables_with_bootstrap": len(bootstrap_quality),
            "invalid_count_before_repair": repaired_total,
            "variables_below_95_percent_finite": len(low_finite),
        }
    )
    if repaired_total:
        audit["bootstrap_quality"]["warnings"].append(
            f"{repaired_total} bootstrap samples were repaired before downstream propagation"
        )
    if low_finite:
        audit["bootstrap_quality"]["warnings"].append("Some variables have less than 95% finite bootstrap samples")

    audit["derived_units"]["checks"].append(
        {
            "a": "a = r0_phys / (r0/a), default r0_phys=0.5 fm",
            "length": "L0 * a",
            "volume_r0": "L0*L1*L2*L3 * a^4",
            "epsilon_bar": "epsilon1 / beta * (r0/a)^2",
        }
    )
    if "derived_summary" not in cache_payload:
        audit["derived_units"]["warnings"].append("No derived summary provided; unit audit is formula-only")

    warning_count = 0
    for section in audit.values():
        if isinstance(section, dict) and isinstance(section.get("warnings"), list):
            warning_count += len(section["warnings"])
    audit["summary"]["warning_count"] = warning_count
    audit["summary"]["v_r_candidate_count"] = len(v_r_records)
    return audit


def load_or_build_scan_cache(
    run_dirs: list[str],
    *,
    cache_path: Path | None = None,
    cache_root: Path | None = None,
    thermalization_steps: int | None = None,
    thermalization_steps_by_run: dict[str, int] | None = None,
    n_bootstrap: int = run_evaluation.DEFAULT_N_BOOTSTRAP,
    block_size: int = run_evaluation.DEFAULT_BLOCK_SIZE,
    target_force: float = run_evaluation.DEFAULT_SOMMER_TARGET,
    r0_t_min: int = run_evaluation.DEFAULT_V_R_T_MIN,
    r0_t_max: int | None = run_evaluation.DEFAULT_V_R_T_MAX,
    load_workers: int = 1,
    calc_workers: int | None = None,
    selected_v_windows: dict[int | str, tuple[int, int | None]] | None = None,
    force_rebuild: bool = False,
) -> dict[str, Any]:
    resolved_run_dirs = _canonical_run_dirs(run_dirs)
    cuts_by_run = resolve_thermalization_steps_by_run(
        resolved_run_dirs,
        thermalization_steps=thermalization_steps,
        thermalization_steps_by_run=thermalization_steps_by_run,
    )
    analysis_options = build_analysis_options(
        thermalization_steps_by_run=cuts_by_run,
        n_bootstrap=n_bootstrap,
        block_size=block_size,
        target_force=target_force,
        r0_t_min=r0_t_min,
        r0_t_max=r0_t_max,
        load_workers=load_workers,
        calc_workers=calc_workers,
        selected_v_windows=selected_v_windows,
    )
    current_input_hash = input_hash_for_runs(resolved_run_dirs)
    expected = expected_cache_metadata(
        resolved_run_dirs,
        analysis_options,
        input_hash=current_input_hash,
    )
    resolved_cache_path = (
        Path(cache_path)
        if cache_path is not None
        else default_scan_cache_path(
            resolved_run_dirs,
            analysis_options,
            cache_root=cache_root,
            input_hash=current_input_hash,
        )
    )

    existing = load_json(resolved_cache_path, default=None)
    stale_reasons = explain_scan_cache_staleness(existing, expected)
    if existing is not None and not force_rebuild and not stale_reasons:
        existing["_cache_status"] = {
            "path": str(resolved_cache_path),
            "loaded_from_cache": True,
            "stale_reasons": [],
        }
        return existing

    combined_w_temp, aggregation = run_evaluation._load_combined_w_temp(
        resolved_run_dirs,
        load_workers=load_workers,
        thermalization_steps_by_run=cuts_by_run,
    )
    if combined_w_temp is None:
        raise RuntimeError("No combined W_temp data could be loaded from the selected runs.")

    calc = Calculator(
        combined_w_temp,
        n_bootstrap=int(n_bootstrap),
        step_size=max(1, int(block_size)),
    )
    unique_rs = [int(r) for r in calc.get_unique_Rs()]
    unique_ts = [int(t) for t in calc.get_unique_Ts()]
    available_pairs = [(int(r), int(t)) for r, t in calc.get_available_pairs()]
    calc.prime_w_rt_cache(available_pairs)

    windows = enumerate_t_windows(unique_ts)
    wrt_scan = build_wrt_scan(calc, unique_rs, unique_ts)
    effective_mass_scan = build_effective_mass_scan(calc, unique_rs, unique_ts)
    v_r_scan = build_v_r_scan(calc, unique_rs, windows)

    selected_v = _selected_v_results(
        calc,
        unique_rs,
        t_min=int(r0_t_min),
        t_max=r0_t_max,
        selected_v_windows=selected_v_windows,
    )
    r0_scan: list[dict[str, Any]] = []
    cornell_curves: list[dict[str, Any]] = []
    if len(selected_v) >= 3:
        r0_scan, cornell_curves = build_locked_r0_scan(
            selected_v,
            target_force=float(target_force),
            n_bootstrap=int(n_bootstrap),
        )
        for row in r0_scan:
            row["t_min"] = int(r0_t_min)
            row["t_max"] = None if r0_t_max is None else int(r0_t_max)
        for row in cornell_curves:
            row["t_min"] = int(r0_t_min)
            row["t_max"] = None if r0_t_max is None else int(r0_t_max)

    payload: dict[str, Any] = {
        **expected,
        "created_at": _utc_now(),
        "aggregation": aggregation,
        "available_R": unique_rs,
        "available_T": unique_ts,
        "available_pairs": [[r, t] for r, t in available_pairs],
        "windows": [
            {"t_min": int(t_min), "t_max": None if t_max is None else int(t_max)}
            for t_min, t_max in windows
        ],
        "wrt_scan": wrt_scan,
        "effective_mass_scan": effective_mass_scan,
        "v_r_scan": v_r_scan,
        "r0_scan": r0_scan,
        "cornell_curves": cornell_curves,
        "bootstrap_quality": _collect_bootstrap_quality(calc),
    }
    payload["audit"] = build_physics_audit(payload)
    save_json(resolved_cache_path, payload)
    payload["_cache_status"] = {
        "path": str(resolved_cache_path),
        "loaded_from_cache": False,
        "stale_reasons": [] if existing is None else stale_reasons,
    }
    return payload


def _get_go():
    try:
        import plotly.graph_objects as go
    except ImportError as exc:
        raise RuntimeError("Plotly is required for notebook analysis figures.") from exc
    return go


def make_wrt_figure(wrt_records: list[dict[str, Any]], R: int | None = None):
    go = _get_go()
    rows = [
        row for row in wrt_records
        if R is None or int(row["R"]) == int(R)
    ]
    fig = go.Figure()
    for r_value in sorted({int(row["R"]) for row in rows}):
        r_rows = sorted([row for row in rows if int(row["R"]) == r_value], key=lambda row: row["T"])
        fig.add_trace(
            go.Scatter(
                x=[int(row["T"]) for row in r_rows],
                y=[float(row["value"]) for row in r_rows],
                mode="lines+markers",
                name=f"R={r_value}",
                error_y={
                    "type": "data",
                    "array": [0.0 if row.get("err") is None else float(row["err"]) for row in r_rows],
                    "visible": any(row.get("err") is not None for row in r_rows),
                },
            )
        )
    fig.update_layout(
        title="Wilson loops W(R,T)",
        xaxis_title="T",
        yaxis_title="W(R,T)",
        template="plotly_white",
    )
    return fig


def make_effective_mass_figure(effective_mass_records: list[dict[str, Any]], R: int | None = None):
    go = _get_go()
    rows = [
        row for row in effective_mass_records
        if R is None or int(row["R"]) == int(R)
    ]
    fig = go.Figure()
    for r_value in sorted({int(row["R"]) for row in rows}):
        r_rows = sorted([row for row in rows if int(row["R"]) == r_value], key=lambda row: row["T"])
        fig.add_trace(
            go.Scatter(
                x=[float(row["t_mid"]) for row in r_rows],
                y=[float(row["value"]) for row in r_rows],
                mode="lines+markers",
                name=f"R={r_value}",
                error_y={
                    "type": "data",
                    "array": [0.0 if row.get("err") is None else float(row["err"]) for row in r_rows],
                    "visible": any(row.get("err") is not None for row in r_rows),
                },
            )
        )
    fig.update_layout(
        title="Effective mass plateau",
        xaxis_title="T + 1/2",
        yaxis_title="m_eff",
        template="plotly_white",
    )
    return fig


def make_vr_window_figure(
    v_r_records: list[dict[str, Any]],
    *,
    t_min: int | None = None,
    t_max: int | None = None,
):
    go = _get_go()
    if t_min is None:
        finite_rows = [row for row in v_r_records if row.get("value") is not None]
        if finite_rows:
            t_min = int(finite_rows[0]["t_min"])
            t_max = finite_rows[0]["t_max"]
    rows = [
        row for row in v_r_records
        if t_min is None
        or (int(row["t_min"]) == int(t_min) and row.get("t_max") == t_max)
    ]
    rows = sorted(rows, key=lambda row: row["R"])
    fig = go.Figure(
        data=[
            go.Scatter(
                x=[int(row["R"]) for row in rows],
                y=[float(row["value"]) for row in rows],
                mode="lines+markers",
                name="V(R)",
                error_y={
                    "type": "data",
                    "array": [0.0 if row.get("err") is None else float(row["err"]) for row in rows],
                    "visible": any(row.get("err") is not None for row in rows),
                },
            )
        ]
    )
    fig.update_layout(
        title=f"Static potential V(R), t_min={t_min}, t_max={t_max}",
        xaxis_title="R",
        yaxis_title="V(R)",
        template="plotly_white",
    )
    return fig


def make_r0_audit_figure(r0_records: list[dict[str, Any]], audit: dict[str, Any] | None = None):
    go = _get_go()
    rows = sorted(r0_records, key=lambda row: row["r_min"])
    fig = go.Figure()
    if rows:
        fig.add_trace(
            go.Scatter(
                x=[int(row["r_min"]) for row in rows],
                y=[float(row["r0"]) for row in rows],
                mode="lines+markers",
                name="r0/a",
                error_y={
                    "type": "data",
                    "array": [0.0 if row.get("err") is None else float(row["err"]) for row in rows],
                    "visible": any(row.get("err") is not None for row in rows),
                },
            )
        )
    warning_count = 0 if audit is None else int(audit.get("summary", {}).get("warning_count", 0))
    fig.update_layout(
        title=f"r0 audit by r_min ({warning_count} warning(s))",
        xaxis_title="r_min",
        yaxis_title="r0/a",
        template="plotly_white",
    )
    return fig


def build_notebook_dashboard_state(cache_payload: dict[str, Any]) -> dict[str, Any]:
    state: dict[str, Any] = {
        "cache": cache_payload,
        "available_R": [int(r) for r in cache_payload.get("available_R", [])],
        "available_T": [int(t) for t in cache_payload.get("available_T", [])],
        "wrt_scan": cache_payload.get("wrt_scan", []),
        "effective_mass_scan": cache_payload.get("effective_mass_scan", []),
        "v_r_scan": cache_payload.get("v_r_scan", []),
        "r0_scan": cache_payload.get("r0_scan", []),
        "audit": cache_payload.get("audit", {}),
    }

    try:
        import ipywidgets as widgets
        from IPython.display import clear_output, display
    except ImportError:
        state["widgets_available"] = False
        return state

    r_options = state["available_R"] or [None]
    t_windows = cache_payload.get("windows", [])
    t_options = [
        (f"t_min={row['t_min']}, t_max={row['t_max']}", (row["t_min"], row["t_max"]))
        for row in t_windows
    ] or [("all", (None, None))]

    r_dropdown = widgets.Dropdown(options=r_options, description="R")
    t_dropdown = widgets.Dropdown(options=t_options, description="Window")
    view_dropdown = widgets.Dropdown(
        options=["W(R,T)", "effective mass", "V(R)", "r0 audit"],
        value="effective mass",
        description="View",
    )
    output = widgets.Output()

    def render(*_args):
        with output:
            clear_output(wait=True)
            view = view_dropdown.value
            r_value = r_dropdown.value
            t_min, t_max = t_dropdown.value
            if view == "W(R,T)":
                display(make_wrt_figure(state["wrt_scan"], R=r_value))
            elif view == "effective mass":
                display(make_effective_mass_figure(state["effective_mass_scan"], R=r_value))
            elif view == "V(R)":
                display(make_vr_window_figure(state["v_r_scan"], t_min=t_min, t_max=t_max))
            else:
                display(make_r0_audit_figure(state["r0_scan"], state["audit"]))

    for control in (r_dropdown, t_dropdown, view_dropdown):
        control.observe(render, names="value")
    render()

    state["widgets_available"] = True
    state["controls"] = {
        "R": r_dropdown,
        "window": t_dropdown,
        "view": view_dropdown,
    }
    state["dashboard"] = widgets.VBox([
        widgets.HBox([view_dropdown, r_dropdown, t_dropdown]),
        output,
    ])
    return state
