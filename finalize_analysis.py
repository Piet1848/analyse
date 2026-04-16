#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import os
import sys
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

import numpy as np

import run_evaluation
from calculator import Calculator, fit_r0_from_potential_data
from finalized_analysis_helpers import (
    build_effective_mass_scan,
    build_locked_r0_scan,
    build_v_r_scan,
    build_wrt_scan,
    enumerate_t_windows,
    load_json,
    save_cornell_plot,
    save_effective_mass_plot,
    save_json,
    save_r0_stability_plot,
    summarize_bootstrap,
    window_label,
)
from load_input_yaml import GaugeObservableParams, MetropolisParams, load_params


SCHEMA_VERSION = 1
IGNORED_METRO_FIELDS = {"seed", "nSweep"}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def format_token(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:g}"
    return str(value)


def analysis_hash(run_dirs: list[str], thermalization_steps: int) -> str:
    payload = {
        "run_dirs": sorted(os.path.abspath(path) for path in run_dirs),
        "thermalization_steps": int(thermalization_steps),
    }
    digest = hashlib.md5(repr(payload).encode("utf-8")).hexdigest()
    return digest[:12]


def compare_dataclass_dicts(reference: dict[str, Any], current: dict[str, Any]) -> dict[str, tuple[Any, Any]]:
    diff: dict[str, tuple[Any, Any]] = {}
    all_keys = sorted(set(reference) | set(current))
    for key in all_keys:
        if reference.get(key) != current.get(key):
            diff[key] = (reference.get(key), current.get(key))
    return diff


def validate_run_directories(
    run_dirs: list[str],
) -> tuple[MetropolisParams, GaugeObservableParams, dict[str, Any]]:
    if not run_dirs:
        raise ValueError("At least one run directory is required.")

    normalized = [os.path.abspath(path) for path in run_dirs]
    missing = [path for path in normalized if not os.path.isdir(path)]
    if missing:
        raise ValueError(f"Run directories not found: {missing}")

    loaded: list[tuple[str, MetropolisParams, GaugeObservableParams]] = []
    for path in normalized:
        yaml_path = Path(path) / "input.yaml"
        if not yaml_path.exists():
            raise ValueError(f"Missing input.yaml in {path}")
        metro, gauge = load_params(str(yaml_path))
        loaded.append((path, metro, gauge))

    ref_path, ref_metro, ref_gauge = loaded[0]
    ref_metro_dict = asdict(ref_metro)
    ref_gauge_dict = asdict(ref_gauge)
    ref_metro_cmp = {k: v for k, v in ref_metro_dict.items() if k not in IGNORED_METRO_FIELDS}

    mismatches: list[str] = []
    for path, metro, gauge in loaded[1:]:
        metro_cmp = {k: v for k, v in asdict(metro).items() if k not in IGNORED_METRO_FIELDS}
        metro_diff = compare_dataclass_dicts(ref_metro_cmp, metro_cmp)
        gauge_diff = compare_dataclass_dicts(ref_gauge_dict, asdict(gauge))
        if not metro_diff and not gauge_diff:
            continue

        parts = [f"Compatibility mismatch for {path} compared to {ref_path}:"]
        for field, (expected, got) in metro_diff.items():
            parts.append(f"  MetropolisParams.{field}: expected {expected!r}, got {got!r}")
        for field, (expected, got) in gauge_diff.items():
            parts.append(f"  GaugeObservableParams.{field}: expected {expected!r}, got {got!r}")
        mismatches.append("\n".join(parts))

    if mismatches:
        raise ValueError("\n".join(mismatches))

    summary = {
        "reference_run": ref_path,
        "ignored_metropolis_fields": sorted(IGNORED_METRO_FIELDS),
        "metropolis_common": ref_metro_cmp,
        "gauge_common": ref_gauge_dict,
    }
    return ref_metro, ref_gauge, summary


class FinalizedAnalysisRunner:
    def __init__(
        self,
        run_dirs: list[str],
        output_root: Path,
        plot_mode: str,
        load_workers: int,
        calc_workers: int | None,
        n_bootstrap: int,
        target_force: float,
        input_func: Callable[[str], str] = input,
        print_func: Callable[..., None] = print,
    ):
        if plot_mode != "html":
            raise ValueError("Only --plot-mode html is currently supported.")

        self.run_dirs = [os.path.abspath(path) for path in run_dirs]
        self.output_root = output_root.resolve()
        self.plot_mode = plot_mode
        self.load_workers = int(load_workers)
        self.calc_workers = calc_workers
        self.n_bootstrap = int(n_bootstrap)
        self.target_force = float(target_force)
        self.input = input_func
        self.print = print_func

        self.metro, self.gauge, self.compatibility_summary = validate_run_directories(self.run_dirs)
        self.thermalization_steps = self._prompt_thermalization()
        self.analysis_id = analysis_hash(self.run_dirs, self.thermalization_steps)
        self.analysis_dir = self._analysis_dir_path()

        self.manifest_path = self.analysis_dir / "manifest.json"
        self.input_runs_path = self.analysis_dir / "input_runs.json"
        self.scan_dir = self.analysis_dir / "scan_cache"
        self.plots_dir = self.analysis_dir / "plots"
        self.wilson_dir = self.analysis_dir / "wilsonloop"
        self.wilson_bootstrap_dir = self.wilson_dir / "bootstrap"
        self.r0_dir = self.analysis_dir / "r0"
        self.r0_bootstrap_dir = self.r0_dir / "bootstrap"
        self.derived_dir = self.analysis_dir / "derived"
        self.derived_bootstrap_dir = self.derived_dir / "bootstrap"

        self.manifest = self._load_or_init_manifest()
        self._ensure_output_layout()
        self._write_input_runs()

        self.combined_w_temp = None
        self.aggregation = None
        self.calc: Calculator | None = None
        self.unique_rs: list[int] = []
        self.unique_ts: list[int] = []
        self.windows: list[tuple[int, int | None]] = []

    def _analysis_dir_path(self) -> Path:
        folder_name = (
            f"beta_{format_token(self.metro.beta)}"
            f"__L_{self.metro.L0}x{self.metro.L1}x{self.metro.L2}x{self.metro.L3}"
            f"__eps1_{format_token(self.metro.epsilon1)}"
            f"__eps2_{format_token(self.metro.epsilon2)}"
            f"__nrun_{len(self.run_dirs)}"
            f"__{self.analysis_id}"
        )
        return self.output_root / folder_name

    def _prompt_thermalization(self) -> int:
        default_cut = int(run_evaluation.THERMALIZATION_STEPS)
        self.print("Selected run directories:")
        for path in self.run_dirs:
            self.print(f"  {path}")
        while True:
            response = self.input(
                f"Thermalization cut is {default_cut}. Press Enter to accept, or type a new integer: "
            ).strip()
            if response == "" or response.lower() in {"y", "yes"}:
                return default_cut
            try:
                value = int(response)
            except ValueError:
                self.print(f"Invalid thermalization cut: {response!r}")
                continue
            if value < 0:
                self.print("Thermalization cut must be non-negative.")
                continue
            return value

    def _load_or_init_manifest(self) -> dict[str, Any]:
        manifest = load_json(self.manifest_path, default=None)
        if manifest is None:
            return {
                "schema_version": SCHEMA_VERSION,
                "analysis_id": self.analysis_id,
                "created_at": utc_now(),
                "updated_at": utc_now(),
                "status": {
                    "scan_cache_built": False,
                    "wilsonloop_complete": False,
                    "r0_complete": False,
                    "derived_complete": False,
                },
                "run_evaluation_version": run_evaluation.CALC_VERSION,
                "calculator_helper": "fit_r0_from_potential_data",
                "thermalization_steps": int(self.thermalization_steps),
                "n_bootstrap": int(self.n_bootstrap),
                "target_force": float(self.target_force),
                "plot_mode": self.plot_mode,
                "calc_workers": self.calc_workers,
                "load_workers": self.load_workers,
                "compatibility_summary": self.compatibility_summary,
            }

        if int(manifest.get("thermalization_steps")) != int(self.thermalization_steps):
            raise ValueError(
                f"Existing analysis uses thermalization_steps={manifest.get('thermalization_steps')} "
                f"but this run requested {self.thermalization_steps}."
            )
        if int(manifest.get("n_bootstrap")) != int(self.n_bootstrap):
            raise ValueError(
                f"Existing analysis uses n_bootstrap={manifest.get('n_bootstrap')} "
                f"but this run requested {self.n_bootstrap}."
            )
        if not np.isclose(float(manifest.get("target_force")), float(self.target_force)):
            raise ValueError(
                f"Existing analysis uses target_force={manifest.get('target_force')} "
                f"but this run requested {self.target_force}."
            )
        manifest["updated_at"] = utc_now()
        return manifest

    def _ensure_output_layout(self) -> None:
        for path in [
            self.analysis_dir,
            self.scan_dir,
            self.plots_dir,
            self.wilson_dir,
            self.wilson_bootstrap_dir,
            self.r0_dir,
            self.r0_bootstrap_dir,
            self.derived_dir,
            self.derived_bootstrap_dir,
        ]:
            path.mkdir(parents=True, exist_ok=True)
        self._save_manifest()

    def _write_input_runs(self) -> None:
        payload = {
            "run_dirs": self.run_dirs,
            "compatibility_summary": self.compatibility_summary,
        }
        save_json(self.input_runs_path, payload)

    def _save_manifest(self) -> None:
        self.manifest["updated_at"] = utc_now()
        save_json(self.manifest_path, self.manifest)

    def _load_data(self) -> None:
        if self.combined_w_temp is not None and self.calc is not None:
            return

        combined_w_temp, aggregation = run_evaluation._load_combined_w_temp(
            self.run_dirs,
            verbose=True,
            prefix="[finalize_analysis]",
            load_workers=self.load_workers,
            thermalization_steps=self.thermalization_steps,
        )
        if combined_w_temp is None:
            raise RuntimeError("No combined W_temp data could be loaded from the selected runs.")

        block_size = max(1, run_evaluation.DEFAULT_BLOCK_SIZE)
        self.combined_w_temp = combined_w_temp
        self.aggregation = aggregation
        self.calc = Calculator(
            combined_w_temp,
            n_bootstrap=self.n_bootstrap,
            step_size=block_size,
        )
        self.unique_rs = [int(r) for r in self.calc.get_unique_Rs()]
        self.unique_ts = [int(t) for t in self.calc.get_unique_Ts()]
        self.windows = enumerate_t_windows(self.unique_ts)

        self.manifest["aggregation"] = aggregation
        self.manifest["block_size"] = block_size
        self.manifest["available_R"] = self.unique_rs
        self.manifest["available_T"] = self.unique_ts
        self._save_manifest()

    def _build_scan_cache(self) -> None:
        self._load_data()
        assert self.calc is not None

        wrt_path = self.scan_dir / "wrt_scan.json"
        eff_path = self.scan_dir / "effective_mass_scan.json"
        vr_path = self.scan_dir / "v_r_scan_candidates.json"

        if not wrt_path.exists():
            save_json(wrt_path, build_wrt_scan(self.calc, self.unique_rs, self.unique_ts))
        if not eff_path.exists():
            save_json(eff_path, build_effective_mass_scan(self.calc, self.unique_rs, self.unique_ts))
        if not vr_path.exists():
            save_json(vr_path, build_v_r_scan(self.calc, self.unique_rs, self.windows))

        self.manifest["status"]["scan_cache_built"] = True
        self._save_manifest()

    def _load_v_choices(self) -> dict[str, Any]:
        return load_json(self.wilson_dir / "choices.json", default={"choices": {}})

    def _load_v_summary(self) -> dict[str, Any]:
        return load_json(self.wilson_dir / "V_R_summary.json", default={"results": {}})

    def _save_v_state(self, choices_payload: dict[str, Any], summary_payload: dict[str, Any]) -> None:
        save_json(self.wilson_dir / "choices.json", choices_payload)
        save_json(self.wilson_dir / "V_R_summary.json", summary_payload)

    def _candidate_windows_for_r(self, r_value: int) -> list[dict[str, Any]]:
        v_scan = load_json(self.scan_dir / "v_r_scan_candidates.json", default=[])
        return [
            row
            for row in v_scan
            if int(row["R"]) == int(r_value)
        ]

    def _prompt_window_choice(self, r_value: int, candidates: list[dict[str, Any]]) -> tuple[int, int | None]:
        if not candidates:
            raise RuntimeError(f"No finite V(R) candidate windows available for R={r_value}.")

        self.print(f"\nR = {r_value}")
        self.print("Candidate windows:")
        for row in candidates:
            self.print(
                f"  {window_label(int(row['t_min']), row['t_max'])}: "
                f"V(R)={row['value']:.6g}"
                + (f" +/- {row['err']:.3g}" if row.get("err") is not None else "")
            )
        while True:
            response = self.input(
                "Choose t_min and optional t_max as `t_min` or `t_min,t_max` (use `none` for no t_max): "
            ).strip()
            if not response:
                self.print("A t_min selection is required.")
                continue

            try:
                parts = [part.strip().lower() for part in response.split(",")]
                if len(parts) == 1:
                    t_min = int(parts[0])
                    t_max = None
                elif len(parts) == 2:
                    t_min = int(parts[0])
                    t_max = None if parts[1] in {"", "none"} else int(parts[1])
                else:
                    raise ValueError
            except ValueError:
                self.print("Please enter either `t_min` or `t_min,t_max`.")
                continue

            for row in candidates:
                if int(row["t_min"]) == t_min and row["t_max"] == t_max:
                    return t_min, t_max

            self.print(f"Window {window_label(t_min, t_max)} is not available for R={r_value}.")

    def _finalize_wilsonloop(self) -> None:
        self._build_scan_cache()
        assert self.calc is not None

        effective_mass_records = load_json(self.scan_dir / "effective_mass_scan.json", default=[])
        v_choices = self._load_v_choices()
        v_summary = self._load_v_summary()

        for r_value in self.unique_rs:
            key = str(int(r_value))
            bootstrap_path = self.wilson_bootstrap_dir / f"V_R_R_{key}.npy"
            if key in v_choices["choices"] and bootstrap_path.exists():
                continue

            candidates = self._candidate_windows_for_r(r_value)
            plot_path = self.plots_dir / f"effective_mass_R_{key}.html"
            save_effective_mass_plot(plot_path, effective_mass_records, candidates, r_value)
            self.print(f"\nSaved effective-mass plot for R={r_value}: {plot_path}")
            t_min, t_max = self._prompt_window_choice(r_value, candidates)

            var = self.calc.get_variable("V_R", R=int(r_value), t_min=int(t_min), t_max=t_max)
            value = var.get()
            if value is None or not np.isfinite(value):
                raise RuntimeError(f"Selected window produced no finite V(R) result for R={r_value}.")

            bootstrap_samples = np.asarray(var.bootstrap(), dtype=float)
            np.save(bootstrap_path, bootstrap_samples)

            choice_payload = {
                "R": int(r_value),
                "t_min": int(t_min),
                "t_max": t_max,
                "plot_path": str(plot_path),
                "saved_at": utc_now(),
            }
            summary_payload = {
                "value": float(value),
                "err": float(var.err()) if var.err() is not None else None,
                "fit_C": (
                    float(var.parameters.get("fit_C"))
                    if var.parameters.get("fit_C") is not None and np.isfinite(var.parameters.get("fit_C"))
                    else None
                ),
                "bootstrap_path": str(bootstrap_path),
            }
            v_choices["choices"][key] = choice_payload
            v_summary["results"][key] = summary_payload
            self._save_v_state(v_choices, v_summary)

        self.manifest["status"]["wilsonloop_complete"] = True
        self._save_manifest()

    def _load_selected_v_results(self) -> dict[int, dict[str, Any]]:
        choices_payload = self._load_v_choices()
        summary_payload = self._load_v_summary()
        results: dict[int, dict[str, Any]] = {}
        for key, choice in choices_payload.get("choices", {}).items():
            summary = summary_payload.get("results", {}).get(key)
            if summary is None:
                continue
            bootstrap_path = Path(summary["bootstrap_path"])
            if not bootstrap_path.is_absolute():
                bootstrap_path = (self.analysis_dir / bootstrap_path).resolve()
            if not bootstrap_path.exists():
                bootstrap_path = self.wilson_bootstrap_dir / f"V_R_R_{key}.npy"
            if not bootstrap_path.exists():
                continue
            results[int(key)] = {
                "R": int(key),
                "t_min": int(choice["t_min"]),
                "t_max": choice["t_max"],
                "value": float(summary["value"]),
                "err": float(summary["err"]) if summary.get("err") is not None else None,
                "fit_C": summary.get("fit_C"),
                "bootstrap_samples": np.load(bootstrap_path),
                "plot_path": choice["plot_path"],
            }
        return results

    def _prompt_r_min_choice(self, r0_records: list[dict[str, Any]]) -> int:
        if not r0_records:
            raise RuntimeError("No candidate r_min values available for the locked V(R) results.")

        self.print("\nCandidate r_min values:")
        for row in sorted(r0_records, key=lambda item: item["r_min"]):
            self.print(
                f"  r_min={row['r_min']}: r0={row['r0']:.6g}"
                + (f" +/- {row['err']:.3g}" if row.get("err") is not None else "")
            )
        valid = {int(row["r_min"]) for row in r0_records}
        while True:
            response = self.input("Choose r_min: ").strip()
            if not response:
                self.print("An r_min selection is required.")
                continue
            try:
                r_min = int(response)
            except ValueError:
                self.print("Please enter an integer r_min.")
                continue
            if r_min not in valid:
                self.print(f"r_min={r_min} is not one of the available candidates: {sorted(valid)}")
                continue
            return r_min

    def _finalize_r0(self) -> None:
        choice_path = self.r0_dir / "choice.json"
        result_path = self.r0_dir / "r0_result.json"
        bootstrap_path = self.r0_bootstrap_dir / "r0.npy"
        if choice_path.exists() and result_path.exists() and bootstrap_path.exists():
            self.manifest["status"]["r0_complete"] = True
            self._save_manifest()
            return

        selected_v = self._load_selected_v_results()
        if len(selected_v) < 3:
            raise RuntimeError("At least three finalized V(R) points are required before fitting r0.")

        r0_records, curve_records = build_locked_r0_scan(
            selected_v,
            target_force=self.target_force,
            n_bootstrap=self.n_bootstrap,
        )
        save_json(
            self.scan_dir / "r0_scan_candidates.json",
            {
                "r0_scan": r0_records,
                "cornell_curves": curve_records,
            },
        )

        stability_path = self.plots_dir / "r0_stability.html"
        cornell_path = self.plots_dir / "cornell_by_rmin.html"
        save_r0_stability_plot(stability_path, r0_records)
        save_cornell_plot(cornell_path, curve_records, r0_records)
        self.print(f"\nSaved r0 stability plot: {stability_path}")
        self.print(f"Saved Cornell fit plot: {cornell_path}")

        r_min = self._prompt_r_min_choice(r0_records)
        fit_rs = [r for r in sorted(selected_v) if r >= int(r_min)]
        rs = np.asarray(fit_rs, dtype=float)
        vs = np.asarray([selected_v[r]["value"] for r in fit_rs], dtype=float)
        errs = np.asarray(
            [selected_v[r]["err"] if selected_v[r]["err"] is not None else 1.0 for r in fit_rs],
            dtype=float,
        )
        boot_matrix = np.stack(
            [np.asarray(selected_v[r]["bootstrap_samples"], dtype=float) for r in fit_rs],
            axis=0,
        )
        fit_result = fit_r0_from_potential_data(
            rs,
            vs,
            errs=errs,
            bootstrap_matrix=boot_matrix,
            target_force=self.target_force,
            n_bootstrap=self.n_bootstrap,
        )

        np.save(bootstrap_path, np.asarray(fit_result["bootstrap_samples"], dtype=float))
        save_json(
            choice_path,
            {
                "r_min": int(r_min),
                "plot_paths": {
                    "stability": str(stability_path),
                    "cornell": str(cornell_path),
                },
                "saved_at": utc_now(),
            },
        )
        save_json(
            self.r0_dir / "cornell_points.json",
            [
                selected_v[r] | {"R": int(r)}
                for r in fit_rs
            ],
        )
        save_json(
            result_path,
            {
                "r0": float(fit_result["r0"]),
                "r0_err": float(fit_result["r0_err"]) if fit_result["r0_err"] is not None else None,
                "target_force": float(self.target_force),
                "cornell_params": fit_result["cornell_params"],
                "fit_rs": [int(r) for r in fit_rs],
                "bootstrap_path": str(bootstrap_path),
            },
        )

        self.manifest["status"]["r0_complete"] = True
        self._save_manifest()

    def _finalize_derived(self) -> None:
        summary_path = self.derived_dir / "summary.json"
        if summary_path.exists():
            self.manifest["status"]["derived_complete"] = True
            self._save_manifest()
            return

        r0_result = load_json(self.r0_dir / "r0_result.json", default=None)
        if r0_result is None:
            raise RuntimeError("r0 result is missing; cannot compute derived quantities.")
        r0_bootstrap = np.load(self.r0_bootstrap_dir / "r0.npy")

        r0_val = float(r0_result["r0"])
        r0_err = r0_result.get("r0_err")
        if not np.isfinite(r0_val) or r0_val <= 0:
            raise RuntimeError("Locked r0 is not finite and positive; cannot compute derived quantities.")

        a_boot = np.full_like(r0_bootstrap, np.nan, dtype=float)
        valid_r0 = np.isfinite(r0_bootstrap) & (r0_bootstrap > 0)
        a_boot[valid_r0] = 0.5 / r0_bootstrap[valid_r0]
        epsilon_boot = np.full_like(r0_bootstrap, np.nan, dtype=float)
        epsilon_boot[valid_r0] = (self.metro.epsilon1 / self.metro.beta) * (r0_bootstrap[valid_r0] ** 2)
        length_boot = np.full_like(r0_bootstrap, np.nan, dtype=float)
        length_boot[valid_r0] = self.metro.L0 * (0.5 / r0_bootstrap[valid_r0])
        volume_boot = np.full_like(r0_bootstrap, np.nan, dtype=float)
        lattice_extent = self.metro.L0 * self.metro.L1 * self.metro.L2 * self.metro.L3
        volume_boot[valid_r0] = lattice_extent * ((0.5 / r0_bootstrap[valid_r0]) ** 4)

        np.save(self.derived_bootstrap_dir / "a.npy", a_boot)
        np.save(self.derived_bootstrap_dir / "epsilon_bar.npy", epsilon_boot)
        np.save(self.derived_bootstrap_dir / "length.npy", length_boot)
        np.save(self.derived_bootstrap_dir / "volume_r0.npy", volume_boot)

        a_mean, a_std = summarize_bootstrap(a_boot)
        eps_mean, eps_std = summarize_bootstrap(epsilon_boot)
        len_mean, len_std = summarize_bootstrap(length_boot)
        vol_mean, vol_std = summarize_bootstrap(volume_boot)

        save_json(
            summary_path,
            {
                "r0": r0_val,
                "r0_err": float(r0_err) if r0_err is not None else None,
                "a": a_mean,
                "a_err": a_std,
                "epsilon_bar": eps_mean,
                "epsilon_bar_err": eps_std,
                "length": len_mean,
                "length_err": len_std,
                "volume_r0": vol_mean,
                "volume_r0_err": vol_std,
                "bootstrap_paths": {
                    "a": str(self.derived_bootstrap_dir / "a.npy"),
                    "epsilon_bar": str(self.derived_bootstrap_dir / "epsilon_bar.npy"),
                    "length": str(self.derived_bootstrap_dir / "length.npy"),
                    "volume_r0": str(self.derived_bootstrap_dir / "volume_r0.npy"),
                },
            },
        )

        self.manifest["status"]["derived_complete"] = True
        self._save_manifest()

    def run(self) -> Path:
        self._finalize_wilsonloop()
        self._finalize_r0()
        self._finalize_derived()
        self.print(f"\nFinalized analysis saved to: {self.analysis_dir}")
        return self.analysis_dir


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Finalize the interactive lattice analysis workflow.")
    parser.add_argument("run_dirs", nargs="+", help="Exact run directories to combine and analyze.")
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("../data/finalized_analysis"),
        help="Root directory for finalized analysis outputs.",
    )
    parser.add_argument(
        "--plot-mode",
        default="html",
        choices=["html"],
        help="Interactive plot mode. Only HTML is currently supported.",
    )
    parser.add_argument("--load-workers", type=int, default=1, help="Parallel workers for loading input files.")
    parser.add_argument("--calc-workers", type=int, default=None, help="Reserved for future parallel calculations.")
    parser.add_argument(
        "--n-bootstrap",
        type=int,
        default=run_evaluation.DEFAULT_N_BOOTSTRAP,
        help="Number of bootstrap replicas for finalized results.",
    )
    parser.add_argument(
        "--target-force",
        type=float,
        default=run_evaluation.DEFAULT_SOMMER_TARGET,
        help="Sommer target force used in the r0 fit.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    runner = FinalizedAnalysisRunner(
        run_dirs=args.run_dirs,
        output_root=args.output_root,
        plot_mode=args.plot_mode,
        load_workers=args.load_workers,
        calc_workers=args.calc_workers,
        n_bootstrap=args.n_bootstrap,
        target_force=args.target_force,
    )
    runner.run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
