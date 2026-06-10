#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Any

import numpy as np

import data_organizer
import finalize_analysis
import run_evaluation
from calculator import Calculator
from finalized_analysis_helpers import (
    _get_plotly,
    _write_figure_html,
    build_bootstrap_block_size_scan,
    recommend_bootstrap_block_size,
    save_bootstrap_block_size_plot,
    save_creutz_diagonal_plot,
    save_gradient_flow_plot,
    save_json,
)


def flow_label(flow_time: float | None) -> str:
    return "unflowed" if flow_time is None else f"{float(flow_time):.12g}"


def filename_token(value: str) -> str:
    return value.replace("-", "m").replace(".", "p")


def save_combined_creutz_flow_diagonal_plot(path: Path, by_flow: dict[str, Any]) -> None:
    go = _get_plotly()
    fig = go.Figure()

    rows_by_flow: dict[float, dict[float, tuple[float, float | None, int]]] = {}
    for _label, payload in by_flow.items():
        if not isinstance(payload, dict):
            continue
        flow_time = payload.get("flow_time") if isinstance(payload, dict) else None
        try:
            plot_flow = 0.0 if flow_time is None else float(flow_time)
        except (TypeError, ValueError):
            continue
        # Unflowed and explicit t_lat=0 are the same curve. If both are present,
        # prefer the explicit t_lat=0 value for duplicate R points.
        source_priority = 0 if flow_time is None else 1
        chi = payload.get("chi", {}) if isinstance(payload, dict) else {}
        chi_err = payload.get("chi_err", {}) if isinstance(payload, dict) else {}
        for key, value in chi.items():
            try:
                r_text, t_text = key.split(",", 1)
                r_val = float(r_text)
                t_val = float(t_text)
                if not np.isclose(r_val, t_val):
                    continue
                err = chi_err.get(key)
                row = (float(value), float(err) if err is not None else None, source_priority)
            except (TypeError, ValueError):
                continue
            current = rows_by_flow.setdefault(plot_flow, {}).get(r_val)
            if current is None or source_priority >= current[2]:
                rows_by_flow[plot_flow][r_val] = row

    flow_values = sorted(flow for flow, rows in rows_by_flow.items() if rows)
    n_flows = len(flow_values)
    if n_flows:
        offset_step = min(0.035, 0.16 / max(n_flows - 1, 1))
    else:
        offset_step = 0.0

    for flow_index, plot_flow in enumerate(flow_values):
        rows = sorted(
            (r_val, value, err)
            for r_val, (value, err, _priority) in rows_by_flow[plot_flow].items()
        )
        if not rows:
            continue
        offset = (flow_index - (n_flows - 1) / 2.0) * offset_step
        fig.add_trace(
            go.Scatter(
                x=[row[0] + offset for row in rows],
                y=[row[1] for row in rows],
                error_y={
                    "type": "data",
                    "array": [row[2] for row in rows],
                    "visible": True,
                },
                customdata=[row[0] for row in rows],
                hovertemplate=(
                    "R=T midpoint=%{customdata:g}<br>"
                    "chi=%{y:g}<extra>t_lat %{fullData.name}</extra>"
                ),
                mode="lines+markers",
                name=f"{plot_flow:g}",
            )
        )

    if not flow_values:
        fig.add_annotation(
            text="No diagonal Creutz-ratio data available",
            x=0.5,
            y=0.5,
            xref="paper",
            yref="paper",
            showarrow=False,
        )

    fig.update_layout(
        title="Creutz ratios on R = T by t_lat",
        template="plotly_white",
        xaxis_title="R = T midpoint",
        yaxis_title="chi(R,T)",
        legend_title="t_lat",
    )
    _write_figure_html(fig, path)


def diagonal_creutz_pairs(available_pairs: list[tuple[int, int]]) -> list[tuple[int, float]]:
    pair_set = {(int(r_val), int(t_val)) for r_val, t_val in available_pairs}
    candidates = sorted({min(r_val, t_val) for r_val, t_val in pair_set})
    result: list[tuple[int, float]] = []
    for base in candidates:
        required = {
            (base, base),
            (base, base + 1),
            (base + 1, base),
            (base + 1, base + 1),
        }
        if required <= pair_set:
            result.append((base, float(base) + 0.5))
    return result


def required_wilson_pairs_for_diagonal(bases: list[int]) -> list[tuple[int, int]]:
    seen: set[tuple[int, int]] = set()
    ordered: list[tuple[int, int]] = []
    for base in bases:
        for pair in (
            (base, base),
            (base, base + 1),
            (base + 1, base),
            (base + 1, base + 1),
        ):
            if pair not in seen:
                ordered.append(pair)
                seen.add(pair)
    return ordered


def build_near_diagonal_thermalization_preview(run_dir: str) -> list[dict[str, Any]]:
    path = Path(run_dir) / "W_temp.out"
    if not path.exists():
        return []

    delimiter: str | None = None
    header_tokens: list[str] | None = None
    header: list[str] | None = None
    steps_by_pair: dict[tuple[int, int], list[float]] = {}
    values_by_pair: dict[tuple[int, int], list[float]] = {}
    pair_order: list[tuple[int, int]] = []

    with path.open("r", encoding="utf-8-sig") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                tokens = stripped[1:].strip().split()
                if tokens and "=" not in stripped and all(data_organizer._looks_like_column_name(token) for token in tokens):
                    header_tokens = tokens
                continue

            if delimiter is None:
                delimiter = "," if "," in stripped else None
            tokens = data_organizer._split_table_line(stripped, delimiter)
            if not tokens:
                continue
            if header is None and all(data_organizer._looks_like_column_name(token) for token in tokens):
                header_tokens = tokens
                continue

            try:
                values = [float(token) for token in tokens]
            except ValueError:
                continue

            if header is None:
                header = data_organizer._canonical_wilson_header(header_tokens, len(values))
                if header is None:
                    return []
            if len(values) != len(header):
                continue
            row = dict(zip(header, values))
            if "t_over_a2" in row:
                continue
            try:
                step = float(row.get("step", row.get("conf_id")))
                r_val = int(row["L"])
                t_val = int(row["T"])
                w_value = float(row["W_temp"])
            except (TypeError, ValueError, KeyError):
                return []
            if not run_evaluation._near_diagonal_wilson_pair_filter(None, r_val, t_val):
                continue

            pair = (r_val, t_val)
            if pair not in values_by_pair:
                pair_order.append(pair)
                steps_by_pair[pair] = []
                values_by_pair[pair] = []
            steps_by_pair[pair].append(step)
            values_by_pair[pair].append(w_value)

    records: list[dict[str, Any]] = []
    run_label = Path(run_dir).name
    for r_val, t_val in pair_order:
        values = np.asarray(values_by_pair[(r_val, t_val)], dtype=float)
        if values.size == 0:
            continue
        running_mean = np.cumsum(values) / np.arange(1, values.size + 1, dtype=float)
        records.append(
            {
                "kind": "run",
                "run_label": run_label,
                "R": int(r_val),
                "T": int(t_val),
                "steps": [float(step) for step in steps_by_pair[(r_val, t_val)]],
                "values": [float(value) for value in values],
                "running_mean": [float(value) for value in running_mean],
            }
        )
    return records


class CreutzFlowAnalysisRunner(finalize_analysis.FinalizedAnalysisRunner):
    def _analysis_dir_path(self) -> Path:
        folder_name = (
            f"creutz_flow__beta_{finalize_analysis.format_token(self.metro.beta)}"
            f"__L_{self.metro.L0}x{self.metro.L1}x{self.metro.L2}x{self.metro.L3}"
            f"__eps1_{finalize_analysis.format_token(self.metro.epsilon1)}"
            f"__nrun_{len(self.run_dirs)}"
            f"__{self.analysis_id}"
        )
        return self.output_root / folder_name

    def _load_or_init_manifest(self) -> dict[str, Any]:
        manifest = super()._load_or_init_manifest()
        manifest["analyser"] = "creutz_flow_diagonal"
        manifest["wilson_loop_filter"] = "near_diagonal_abs_R_minus_T_le_1"
        manifest.setdefault("status", {})["creutz_flow_complete"] = False
        return manifest

    def _thermalization_preview_records_for_run(self, run_dir: str) -> list[dict[str, Any]]:
        run_dir = os.path.abspath(run_dir)
        cached = self.thermalization_preview_records_by_run.get(run_dir)
        if cached is None:
            cached = build_near_diagonal_thermalization_preview(run_dir)
            self.thermalization_preview_records_by_run[run_dir] = cached
        return cached

    def _select_bootstrap_block_size_for_flow_scan(self, file_data) -> int:
        saved_block_size = self._saved_bootstrap_block_size()
        if saved_block_size is not None:
            self.print(f"Reusing bootstrap block size: {saved_block_size}")
            self.manifest["block_size"] = int(saved_block_size)
            self._save_manifest()
            return int(saved_block_size)

        probe_calc = Calculator(file_data, n_bootstrap=1, step_size=1)
        available_flow_times = probe_calc.get_available_flow_times()
        scan_flow_time = None if None in available_flow_times else (available_flow_times[0] if available_flow_times else None)

        self.print("\nBuilding bootstrap block-size diagnostic for representative diagonal W(R,T) points...")
        block_size_scan_values = self._resolved_bootstrap_block_sizes(file_data)
        self.manifest["block_size_scan_values"] = block_size_scan_values
        self._save_manifest()

        scan_records, selected_pairs = build_bootstrap_block_size_scan(
            file_data,
            block_size_scan_values,
            n_bootstrap=self.n_bootstrap,
            flow_time=scan_flow_time,
        )
        tau_hint = self._compute_tau_hint(file_data)
        recommended = recommend_bootstrap_block_size(
            scan_records,
            fallback=run_evaluation.DEFAULT_BLOCK_SIZE,
        )

        scan_path = self.scan_dir / "bootstrap_block_size_scan.json"
        plot_path = self.plots_dir / "bootstrap_block_size.html"
        save_json(
            scan_path,
            {
                "schema_version": finalize_analysis.SCHEMA_VERSION,
                "block_size_scan_values": block_size_scan_values,
                "block_size_scan_total_length": self._bootstrap_scan_total_length(file_data),
                "block_size_scan_flow_time": None if scan_flow_time is None else float(scan_flow_time),
                "selected_pairs": [{"R": int(r_val), "T": int(t_val)} for r_val, t_val in selected_pairs],
                "tau_int": tau_hint,
                "recommended_block_size": int(recommended),
                "records": scan_records,
                "saved_at": finalize_analysis.utc_now(),
            },
        )
        save_bootstrap_block_size_plot(
            plot_path,
            scan_records,
            recommended_block_size=int(recommended),
        )
        opened = self._open_plot(plot_path)
        action = "Opened" if opened else "Saved"
        self.print(f"{action} bootstrap block-size plot: {plot_path}")

        block_size = self._prompt_bootstrap_block_size(int(recommended))
        save_json(
            self.bootstrap_block_size_choice_path,
            {
                "schema_version": finalize_analysis.SCHEMA_VERSION,
                "block_size": int(block_size),
                "recommended_block_size": int(recommended),
                "scan_path": str(scan_path),
                "plot_path": str(plot_path),
                "block_size_scan_flow_time": None if scan_flow_time is None else float(scan_flow_time),
                "tau_int": tau_hint,
                "saved_at": finalize_analysis.utc_now(),
            },
        )
        self.manifest["block_size"] = int(block_size)
        self.manifest["bootstrap_block_size_scan"] = str(scan_path)
        self.manifest["bootstrap_block_size_plot"] = str(plot_path)
        self.manifest["block_size_scan_flow_time"] = None if scan_flow_time is None else float(scan_flow_time)
        self.manifest["tau_int"] = tau_hint
        self._save_manifest()
        return int(block_size)

    def _load_filtered_data(self) -> int:
        if self.combined_w_temp is not None:
            return int(self.manifest.get("block_size") or 1)

        combined_w_temp, aggregation = run_evaluation._load_combined_w_temp_filtered(
            self.run_dirs,
            pair_filter=run_evaluation._near_diagonal_wilson_pair_filter,
            verbose=True,
            prefix="[analyze_creutz_flow]",
            load_workers=self.load_workers,
            thermalization_steps_by_run=self.thermalization_steps_by_run,
        )
        combined_gradient_flow, gradient_flow_metadata = run_evaluation._load_combined_gradient_flow_obs(
            self.run_dirs,
            load_workers=self.load_workers,
            thermalization_steps_by_run=self.thermalization_steps_by_run,
        )
        if combined_w_temp is None:
            raise RuntimeError("No near-diagonal W_temp data could be loaded from the selected runs.")

        self.print(f"Combined filtered W_temp data loaded with aggregation: {aggregation}")
        block_size = self._select_bootstrap_block_size_for_flow_scan(combined_w_temp)

        probe_calc = Calculator(combined_w_temp, n_bootstrap=1, step_size=1)
        available_flow_times = probe_calc.get_available_flow_times()
        self.combined_w_temp = combined_w_temp
        self.combined_gradient_flow = combined_gradient_flow
        self.gradient_flow_metadata = gradient_flow_metadata
        self.aggregation = aggregation
        self.manifest["aggregation"] = aggregation
        self.manifest["block_size"] = block_size
        self.manifest["available_wilson_flow_times"] = [
            None if value is None else float(value)
            for value in available_flow_times
        ]
        self.manifest["gradient_flow_metadata"] = gradient_flow_metadata
        self._save_manifest()
        return int(block_size)

    def _finalize_gradient_flow(self) -> None:
        block_size = self._load_filtered_data()
        summary = run_evaluation.summarize_gradient_flow_obs(
            self.combined_gradient_flow,
            t0_target=self.gradient_flow_metadata.get("t0_target"),
            block_size=block_size,
            n_bootstrap=self.n_bootstrap,
            include_bootstrap=True,
        )
        bootstrap_samples = summary.pop("bootstrap_samples", {}) if summary else {}
        summary_path = self.gradient_flow_dir / "summary.json"
        plot_path = self.plots_dir / "gradient_flow_t2E.html"

        if summary:
            for name, samples in bootstrap_samples.items():
                np.save(self.gradient_flow_bootstrap_dir / f"{name}.npy", np.asarray(samples, dtype=float))
            summary["bootstrap_paths"] = {
                name: str(self.gradient_flow_bootstrap_dir / f"{name}.npy")
                for name in bootstrap_samples
            }
            save_json(summary_path, summary)
            save_gradient_flow_plot(plot_path, summary)
            self.manifest["gradient_flow_summary"] = str(summary_path)
            self.manifest["gradient_flow_plot"] = str(plot_path)
        else:
            save_json(summary_path, {"status": "no gradient_flow_obs.dat data available"})
            self.manifest["gradient_flow_summary"] = str(summary_path)

        self.manifest["status"]["gradient_flow_complete"] = True
        self._save_manifest()

    def _finalize_creutz_flow(self) -> None:
        block_size = self._load_filtered_data()
        assert self.combined_w_temp is not None

        probe_calc = Calculator(self.combined_w_temp, n_bootstrap=1, step_size=1)
        available_flow_times = probe_calc.get_available_flow_times()
        by_flow: dict[str, Any] = {}
        plot_paths: dict[str, str] = {}
        bootstrap_root = self.creutz_bootstrap_dir / "diagonal_by_flow_time"
        bootstrap_root.mkdir(parents=True, exist_ok=True)

        for flow_time in available_flow_times:
            label = flow_label(flow_time)
            token = filename_token(label)
            calc = Calculator(self.combined_w_temp, n_bootstrap=self.n_bootstrap, step_size=block_size)
            available_pairs = calc.get_available_pairs(flow_time=flow_time)
            diagonal_pairs = diagonal_creutz_pairs(available_pairs)
            bases = [base for base, _ in diagonal_pairs]
            required_pairs = required_wilson_pairs_for_diagonal(bases)
            if required_pairs:
                calc.prime_w_rt_cache(pairs=required_pairs, flow_time=flow_time)

            chi: dict[str, float] = {}
            chi_err: dict[str, float] = {}
            chi_boot_paths: dict[str, str] = {}
            flow_bootstrap_dir = bootstrap_root / token
            flow_bootstrap_dir.mkdir(parents=True, exist_ok=True)

            for _base, midpoint in diagonal_pairs:
                try:
                    var = calc.get_variable("chi", R=midpoint, T=midpoint, flow_time=flow_time)
                except Exception:
                    continue
                value = var.get()
                if value is None or not np.isfinite(value):
                    continue
                key = f"{midpoint:g},{midpoint:g}"
                chi[key] = float(value)
                err = var.err()
                if err is not None and np.isfinite(err):
                    chi_err[key] = float(err)
                boot = var.bootstrap()
                if boot is not None:
                    path = flow_bootstrap_dir / f"chi_R_eq_T_{filename_token(f'{midpoint:g}')}.npy"
                    np.save(path, np.asarray(boot, dtype=float))
                    chi_boot_paths[key] = str(path)

            flow_payload = {
                "status": "ok" if chi else "not enough adjacent diagonal Wilson-loop data",
                "flow_time": None if flow_time is None else float(flow_time),
                "available_pairs_loaded": len(available_pairs),
                "required_pairs_used": [{"R": int(r), "T": int(t)} for r, t in required_pairs],
                "chi": chi,
                "chi_err": chi_err,
                "chi_bootstrap_paths": chi_boot_paths,
            }
            by_flow[label] = flow_payload

            plot_path = self.plots_dir / f"creutz_ratios_R_eq_T__flow_{token}.html"
            save_creutz_diagonal_plot(plot_path, flow_payload)
            plot_paths[label] = str(plot_path)

        combined_plot_path = self.plots_dir / "creutz_ratios_R_eq_T__all_flow_times.html"
        save_combined_creutz_flow_diagonal_plot(combined_plot_path, by_flow)

        summary_path = self.creutz_dir / "diagonal_flow_summary.json"
        payload = {
            "status": "ok" if any(flow.get("chi") for flow in by_flow.values()) else "no diagonal Creutz-ratio data available",
            "available_wilson_flow_times": [
                None if value is None else float(value)
                for value in available_flow_times
            ],
            "n_bootstrap": int(self.n_bootstrap),
            "block_size": int(block_size),
            "thermalization_steps": self.thermalization_steps,
            "thermalization_steps_by_run": self.thermalization_steps_by_run,
            "wilson_loop_filter": "near_diagonal_abs_R_minus_T_le_1",
            "by_flow_time": by_flow,
            "plot_paths": plot_paths,
            "combined_plot_path": str(combined_plot_path),
        }
        save_json(summary_path, payload)
        self.manifest["creutz_diagonal_flow_summary"] = str(summary_path)
        self.manifest["creutz_diagonal_flow_plots"] = plot_paths
        self.manifest["creutz_diagonal_flow_combined_plot"] = str(combined_plot_path)
        self.manifest["status"]["creutz_flow_complete"] = True
        self._save_manifest()

    def run(self) -> Path:
        self._finalize_gradient_flow()
        self._finalize_creutz_flow()
        self.print(f"\nCreutz-flow analysis saved to: {self.analysis_dir}")
        self.print(f"Saved diagonal Creutz summary: {self.creutz_dir / 'diagonal_flow_summary.json'}")
        self.print(f"Saved gradient-flow t2E summary: {self.gradient_flow_dir / 'summary.json'}")
        return self.analysis_dir


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Analyze diagonal Creutz ratios for every Wilson-loop t_lat value and gradient-flow t2E."
    )
    parser.add_argument("run_dirs", nargs="*", help="Exact run directories, or roots containing run directories.")
    parser.add_argument("--run-root", nargs="+", help="Discover run directories beneath one or more roots.")
    parser.add_argument("--filter", nargs="+", metavar="NAME=VALUE", help="YAML parameter filters, e.g. beta=2.4 eps1=0.0 L0=24.")
    parser.add_argument("--exclude", nargs="+", metavar="TEXT_OR_GLOB", help="Exclude run directories whose path contains this text or matches this glob.")
    parser.add_argument("--group-by", nargs="+", help="Optional field names used to split discovered runs.")
    parser.add_argument("--require-fixed-dt", action="store_true", help="Reject or split combined runs whose GradientFlowParams.dt values differ.")
    parser.add_argument("--min-group-size", type=int, default=1, help="Skip discovered groups smaller than this.")
    parser.add_argument("--list-groups", "--list-group", action="store_true", help="Print discovered groups and exit.")
    parser.add_argument("--output-root", type=Path, default=Path("../data/creutz_flow_analysis"), help="Output root.")
    parser.add_argument("--plot-mode", default="html", choices=["html"], help="Only HTML is currently supported.")
    parser.add_argument("--load-workers", type=int, default=1, help="Parallel workers for loading input files.")
    parser.add_argument("--n-bootstrap", type=int, default=run_evaluation.DEFAULT_N_BOOTSTRAP, help="Number of bootstrap replicas.")
    parser.add_argument("--block-sizes", help="Comma-separated bootstrap block sizes.")
    parser.add_argument("--min", "--min-block-size", dest="min_block_size", type=int, default=finalize_analysis.DEFAULT_BLOCK_SIZE_SCAN_MIN)
    parser.add_argument("--max", "--max-block-size", dest="max_block_size", type=int, default=finalize_analysis.DEFAULT_BLOCK_SIZE_SCAN_MAX)
    parser.add_argument("--block-step", type=int, default=finalize_analysis.DEFAULT_BLOCK_SIZE_SCAN_STEP)
    parser.add_argument("--no-open-plots", action="store_false", dest="open_plots", help="Do not auto-open HTML diagnostics.")
    parser.set_defaults(open_plots=True)
    return parser


def grouped_runs_from_args(args: argparse.Namespace) -> list[tuple[tuple[tuple[str, Any], ...], list[str]]]:
    if args.run_root and args.run_dirs:
        raise ValueError("Use either positional run_dirs or --run-root roots, not both.")
    if not args.run_root and not args.run_dirs:
        raise ValueError("Provide run directories or use --run-root.")
    if args.min_group_size < 1:
        raise ValueError("--min-group-size must be at least 1.")

    filter_criteria = finalize_analysis.parse_filter_tokens(args.filter)

    if args.run_root:
        run_dirs = finalize_analysis.discover_run_directories_from_roots(args.run_root)
        run_dirs = finalize_analysis.filter_run_directories(run_dirs, filter_criteria)
        run_dirs = finalize_analysis.exclude_run_directories(run_dirs, args.exclude)
        expanded_roots = True
    else:
        run_dirs, expanded_roots = finalize_analysis.resolve_selected_run_directories(args.run_dirs)
        run_dirs = finalize_analysis.filter_run_directories(run_dirs, filter_criteria)
        run_dirs = finalize_analysis.exclude_run_directories(run_dirs, args.exclude)

    if not run_dirs:
        raise ValueError("No selected run directories remain.")

    if args.run_root or expanded_roots:
        grouped = finalize_analysis.group_run_directories(
            run_dirs,
            group_by=args.group_by,
            require_fixed_dt=args.require_fixed_dt,
        )
    else:
        grouped = [(tuple(), run_dirs)]

    grouped = [
        (signature, paths)
        for signature, paths in grouped
        if len(paths) >= args.min_group_size
    ]
    if not grouped:
        raise ValueError(f"No discovered groups have at least {args.min_group_size} matching run(s).")
    return grouped


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        block_size_scan_values = finalize_analysis.parse_bootstrap_block_sizes(args)
        grouped_runs = grouped_runs_from_args(args)
    except ValueError as exc:
        parser.error(str(exc))

    if args.list_groups:
        for index, (signature, run_dirs) in enumerate(grouped_runs, start=1):
            header = finalize_analysis.describe_group_signature(signature) if signature else "manual selection"
            print(f"[group {index}] {header}")
            for run_dir in run_dirs:
                print(f"  {run_dir}")
        return 0

    for index, (signature, run_dirs) in enumerate(grouped_runs, start=1):
        if len(grouped_runs) > 1:
            header = finalize_analysis.describe_group_signature(signature)
            print(f"\n=== Creutz-flow analysis group {index}/{len(grouped_runs)}: {header} ===")
        runner = CreutzFlowAnalysisRunner(
            run_dirs=run_dirs,
            output_root=args.output_root,
            plot_mode=args.plot_mode,
            load_workers=args.load_workers,
            calc_workers=None,
            n_bootstrap=args.n_bootstrap,
            target_force=run_evaluation.DEFAULT_SOMMER_TARGET,
            block_size_scan_values=block_size_scan_values,
            block_size_scan_min=args.min_block_size,
            block_size_scan_max=args.max_block_size,
            block_size_scan_step=args.block_step,
            wilson_flow_time=None,
            require_fixed_dt=args.require_fixed_dt,
            open_plots=args.open_plots,
        )
        runner.run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
