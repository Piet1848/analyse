from __future__ import annotations

import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import yaml

import finalize_analysis
import finalized_analysis_helpers as helpers
import run_evaluation
import search_data
import data_organizer
from calculator import Calculator, creutz_chi_from_wilson, fit_r0_from_potential_data
from load_input_yaml import load_params


def write_input_yaml(
    run_dir: Path,
    *,
    beta: float = 2.4,
    epsilon1: float = 0.2,
    epsilon2: float = 0.0,
    seed: int = 1,
    n_sweep: int = 100,
    gradient_dt: float = 0.05,
    gradient_t_values: list[float] | None = None,
    w_temp_pairs=None,
) -> None:
    payload = {
        "MetropolisParams": {
            "Ndims": 4,
            "Nd": 4,
            "Nc": 2,
            "L0": 4,
            "L1": 4,
            "L2": 4,
            "L3": 4,
            "nHits": 10,
            "nSweep": n_sweep,
            "seed": seed,
            "beta": beta,
            "delta": 0.1,
            "epsilon1": epsilon1,
            "epsilon2": epsilon2,
        },
        "GaugeObservableParams": {
            "measurement_interval": 1,
            "measure_plaquette": False,
            "measure_wilson_loop_temporal": True,
            "measure_wilson_loop_mu_nu": False,
            "measure_retrace_U": False,
            "W_temp_L_T_pairs": [[2, "1:4"], [3, "1:4"], [4, "1:4"]] if w_temp_pairs is None else w_temp_pairs,
            "W_mu_nu_pairs": [],
            "W_Lmu_Lnu_pairs": [],
            "plaquette_filename": "plaquette.out",
            "W_temp_filename": "W_temp.out",
            "W_mu_nu_filename": "W_mu_nu.out",
            "RetraceU_filename": "RetraceU.out",
            "write_to_file": True,
        },
        "GradientFlowParams": {
            "enabled": True,
            "integrator": "rk3",
            "dt": gradient_dt,
            "t_values": [0.0, 0.25] if gradient_t_values is None else gradient_t_values,
            "measure_energy_clover": True,
            "measure_wilson_loop_temporal": False,
            "measure_wilson_loop_mu_nu": False,
            "extract_t0": False,
            "t0_target": 0.1,
            "obs_filename": "gradient_flow_obs.dat",
            "W_temp_filename": "gradient_flow_wtemp.dat",
            "W_mu_nu_filename": "gradient_flow_w_mu_nu.dat",
            "t0_filename": "gradient_flow_t0.dat",
        },
    }
    with (run_dir / "input.yaml").open("w", encoding="utf-8") as handle:
        yaml.safe_dump(payload, handle, sort_keys=False)


def write_w_temp(
    run_dir: Path,
    *,
    run_offset: float = 0.0,
    r_values: list[int] | None = None,
    t_values: list[int] | None = None,
) -> None:
    A = 0.5
    sigma = 0.2
    B = 0.1
    steps = [0, 1000, 2000, 3000]
    r_values = [2, 3, 4] if r_values is None else r_values
    t_values = [1, 2, 3, 4] if t_values is None else t_values
    rows = ["step,L,T,W_temp"]
    for cfg_index, step in enumerate(steps):
        for r_val in r_values:
            potential = A + sigma * r_val - B / r_val
            for t_val in t_values:
                modulation = 1.0 + 0.01 * cfg_index + run_offset
                value = modulation * np.exp(-potential * t_val)
                rows.append(f"{step},{r_val},{t_val},{value:.10f}")
    (run_dir / "W_temp.out").write_text("\n".join(rows) + "\n", encoding="utf-8")


def write_gradient_flow_obs(
    run_dir: Path,
    rows: list[tuple[int, float, float, float]],
) -> None:
    lines = ["step t_over_a2 Ehat_clover t2E_clover"]
    lines.extend(
        f"{step} {t_over_a2:.12g} {ehat:.12g} {t2e:.12g}"
        for step, t_over_a2, ehat, t2e in rows
    )
    (run_dir / "gradient_flow_obs.dat").write_text("\n".join(lines) + "\n", encoding="utf-8")


def make_run(
    base_dir: Path,
    name: str,
    *,
    beta: float = 2.4,
    epsilon1: float = 0.2,
    epsilon2: float = 0.0,
    seed: int = 1,
    n_sweep: int = 100,
    run_offset: float = 0.0,
    gradient_dt: float = 0.05,
    gradient_t_values: list[float] | None = None,
    w_temp_pairs=None,
    w_temp_r_values: list[int] | None = None,
    w_temp_t_values: list[int] | None = None,
) -> Path:
    run_dir = base_dir / name
    run_dir.mkdir(parents=True, exist_ok=True)
    write_input_yaml(
        run_dir,
        beta=beta,
        epsilon1=epsilon1,
        epsilon2=epsilon2,
        seed=seed,
        n_sweep=n_sweep,
        gradient_dt=gradient_dt,
        gradient_t_values=gradient_t_values,
        w_temp_pairs=w_temp_pairs,
    )
    write_w_temp(
        run_dir,
        run_offset=run_offset,
        r_values=w_temp_r_values,
        t_values=w_temp_t_values,
    )
    return run_dir


def write_stub_plot(path: Path, *args, **kwargs) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("<html><body>stub</body></html>", encoding="utf-8")


def write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


class InputSequence:
    def __init__(self, responses):
        self._responses = iter(responses)

    def __call__(self, prompt: str) -> str:
        try:
            return next(self._responses)
        except StopIteration as exc:
            raise RuntimeError("input exhausted") from exc


class FinalizeAnalysisTests(unittest.TestCase):
    def test_compact_result_summary_formatting_matches_requested_precision(self):
        self.assertEqual(
            finalize_analysis.format_result_value_rounded_to_error(4.38115, 0.0156),
            "4.381 +/- 0.016",
        )
        self.assertEqual(
            finalize_analysis.format_result_value_rounded_to_error(2.73925, 0.00977),
            "2.739 +/- 0.010",
        )
        self.assertEqual(finalize_analysis.format_compact_float(0.0), "0.0")
        self.assertEqual(
            finalize_analysis.format_result_value_rounded_to_error(1.84789, 0.000638),
            "1.8479 +/- 0.0006",
        )

    def test_prompt_window_choice_reuses_previous_on_blank_input(self):
        runner = finalize_analysis.FinalizedAnalysisRunner.__new__(finalize_analysis.FinalizedAnalysisRunner)
        runner.input = InputSequence([""])
        runner.print = lambda *args, **kwargs: None

        choice = runner._prompt_window_choice(
            3,
            [{"R": 3, "t_min": 1, "t_max": None, "value": 0.1}],
            previous_choice=(1, None),
        )

        self.assertEqual(choice, (1, None))

    def test_select_preview_pair_keys_limits_to_representative_subset(self):
        selected = helpers._select_preview_pair_keys(
            [(1, 1), (1, 8), (5, 8), (8, 10), (8, 11), (10, 10)],
        )

        self.assertEqual(selected, [(1, 1), (1, 8), (5, 8), (8, 10), (10, 10)])

    def test_trim_preview_series_keeps_first_ten_percent(self):
        row = {
            "steps": list(range(100)),
            "values": [float(idx) for idx in range(100)],
            "running_mean": [float(idx) / 2.0 for idx in range(100)],
        }

        trimmed = helpers._trim_preview_series(row, fraction=0.10)

        self.assertEqual(trimmed["steps"], list(range(10)))
        self.assertEqual(trimmed["values"], [float(idx) for idx in range(10)])
        self.assertEqual(trimmed["running_mean"], [float(idx) / 2.0 for idx in range(10)])

    def test_prompt_bootstrap_block_size_rejects_invalid_values(self):
        runner = finalize_analysis.FinalizedAnalysisRunner.__new__(finalize_analysis.FinalizedAnalysisRunner)
        runner.input = InputSequence(["abc", "0", "3"])
        messages = []
        runner.print = lambda *args, **kwargs: messages.append(" ".join(str(arg) for arg in args))

        block_size = runner._prompt_bootstrap_block_size(2)

        self.assertEqual(block_size, 3)
        self.assertTrue(any("Invalid bootstrap block size" in message for message in messages))
        self.assertTrue(any("must be at least 1" in message for message in messages))

    def test_prompt_thermalization_star_applies_to_all_runs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_a = make_run(root, "run_a")
            run_b = make_run(root, "run_b", seed=2, n_sweep=200)

            runner = finalize_analysis.FinalizedAnalysisRunner.__new__(finalize_analysis.FinalizedAnalysisRunner)
            runner.run_dirs = [str(run_a.resolve()), str(run_b.resolve())]
            runner.input = InputSequence(["2200*"])
            runner.print = lambda *args, **kwargs: None
            runner._load_saved_thermalization_cuts = lambda: {}
            runner._save_selection_thermalization_preview = lambda run_dir, suggested_cut: None
            runner._save_thermalization_selection_state = lambda cuts_by_run: None

            cuts = runner._prompt_thermalization_by_run()

            self.assertEqual(
                cuts,
                {
                    str(run_a.resolve()): 2200,
                    str(run_b.resolve()): 2200,
                },
            )

    def test_default_bootstrap_scan_extends_to_one_sixteenth_of_total_length(self):
        runner = finalize_analysis.FinalizedAnalysisRunner.__new__(finalize_analysis.FinalizedAnalysisRunner)
        runner.block_size_scan_values = None
        runner.block_size_scan_min = 1
        runner.block_size_scan_max = None
        runner.block_size_scan_step = 1

        file_data = mock.Mock(n_configurations=80_000)

        block_sizes = runner._resolved_bootstrap_block_sizes(file_data)

        self.assertEqual(block_sizes[0], 1)
        self.assertEqual(block_sizes[-1], 5_000)

    def test_bootstrap_scan_includes_dynamic_max_with_coarse_step(self):
        runner = finalize_analysis.FinalizedAnalysisRunner.__new__(finalize_analysis.FinalizedAnalysisRunner)
        runner.block_size_scan_values = None
        runner.block_size_scan_min = 1
        runner.block_size_scan_max = None
        runner.block_size_scan_step = 750

        file_data = mock.Mock(n_configurations=80_000)

        block_sizes = runner._resolved_bootstrap_block_sizes(file_data)

        self.assertEqual(block_sizes[-1], 5_000)

    def test_open_html_plot_returns_false_when_linux_openers_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            html_path = Path(tmp) / "preview.html"
            html_path.write_text("<html></html>", encoding="utf-8")

            failing_process = mock.Mock()
            failing_process.wait.return_value = 1

            with (
                mock.patch.object(finalize_analysis.os, "name", "posix"),
                mock.patch.object(finalize_analysis.sys, "platform", "linux"),
                mock.patch.dict(finalize_analysis.os.environ, {}, clear=True),
                mock.patch.object(
                    finalize_analysis.shutil,
                    "which",
                    side_effect=lambda name: {
                        "xdg-open": "/usr/bin/xdg-open",
                        "gio": "/usr/bin/gio",
                    }.get(name),
                ),
                mock.patch.object(finalize_analysis.subprocess, "Popen", return_value=failing_process) as popen_mock,
            ):
                opened = finalize_analysis.open_html_plot(html_path)

            self.assertFalse(opened)
            self.assertEqual(popen_mock.call_count, 2)
            first_call = popen_mock.call_args_list[0]
            self.assertEqual(first_call.args[0], ["/usr/bin/xdg-open", str(html_path.resolve())])
            self.assertIs(first_call.kwargs["stdin"], finalize_analysis.subprocess.DEVNULL)
            self.assertIs(first_call.kwargs["stdout"], finalize_analysis.subprocess.DEVNULL)
            self.assertIs(first_call.kwargs["stderr"], finalize_analysis.subprocess.DEVNULL)
            self.assertTrue(first_call.kwargs["start_new_session"])

    def test_validate_run_directories_accepts_seed_and_nsweep_differences(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_a = make_run(root, "run_a", seed=1, n_sweep=100)
            run_b = make_run(root, "run_b", seed=2, n_sweep=200, run_offset=0.005)

            metro, gauge, summary = finalize_analysis.validate_run_directories([str(run_a), str(run_b)])

            self.assertAlmostEqual(metro.beta, 2.4)
            self.assertTrue(gauge.measure_wilson_loop_temporal)
            self.assertEqual(summary["ignored_metropolis_fields"], ["nSweep", "seed"])

    def test_validate_run_directories_rejects_non_ignored_mismatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_a = make_run(root, "run_a", beta=2.4)
            run_b = make_run(root, "run_b", beta=2.6)

            with self.assertRaisesRegex(ValueError, r"MetropolisParams\.beta"):
                finalize_analysis.validate_run_directories([str(run_a), str(run_b)])

    def test_validate_run_directories_allows_measurement_grid_mismatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_a = make_run(root, "run_a", w_temp_pairs=[[2, "1:4"], [3, "1:4"]])
            run_b = make_run(root, "run_b", seed=2, w_temp_pairs=[[2, "1:6"], [3, "1:6"], [4, "1:6"]])

            _, _, summary = finalize_analysis.validate_run_directories([str(run_a), str(run_b)])

            self.assertIn("W_temp_L_T_pairs", summary["gauge_varying"])
            self.assertAlmostEqual(summary["metropolis_common"]["beta"], 2.4)

    def test_validate_run_directories_allows_gradient_flow_dt_and_t_values_mismatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_a = make_run(root, "run_a", gradient_dt=0.05, gradient_t_values=[0.0, 0.25])
            run_b = make_run(root, "run_b", seed=2, gradient_dt=0.25, gradient_t_values=[0.0, 0.5])

            _, _, summary = finalize_analysis.validate_run_directories([str(run_a), str(run_b)])

            self.assertEqual(summary["gradient_flow_dt_values"], [0.05, 0.25])
            self.assertIn("dt", summary["gradient_flow_varying"])
            self.assertIn("t_values", summary["gradient_flow_varying"])

    def test_validate_run_directories_can_require_fixed_gradient_flow_dt(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_a = make_run(root, "run_a", gradient_dt=0.05)
            run_b = make_run(root, "run_b", seed=2, gradient_dt=0.25)

            with self.assertRaisesRegex(ValueError, r"GradientFlowParams\.dt"):
                finalize_analysis.validate_run_directories(
                    [str(run_a), str(run_b)],
                    require_fixed_dt=True,
                )

    def test_group_run_directories_combines_gradient_flow_dt_by_default(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_a = make_run(root, "run_a", gradient_dt=0.05)
            run_b = make_run(root, "run_b", seed=2, gradient_dt=0.05)
            run_c = make_run(root, "run_c", seed=3, gradient_dt=0.25)

            groups = finalize_analysis.group_run_directories(
                finalize_analysis.discover_run_directories(str(root)),
            )

            grouped_paths = [paths for _, paths in groups]
            self.assertEqual(
                grouped_paths,
                [sorted([str(run_a.resolve()), str(run_b.resolve()), str(run_c.resolve())])],
            )

    def test_group_run_directories_combines_measurement_grid_by_default(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_a = make_run(root, "run_a", w_temp_pairs=[[2, "1:4"], [3, "1:4"]])
            run_b = make_run(root, "run_b", seed=2, w_temp_pairs=[[2, "1:6"], [3, "1:6"], [4, "1:6"]])

            groups = finalize_analysis.group_run_directories(
                finalize_analysis.discover_run_directories(str(root)),
            )

            grouped_paths = [paths for _, paths in groups]
            self.assertEqual(grouped_paths, [sorted([str(run_a.resolve()), str(run_b.resolve())])])

    def test_group_run_directories_splits_by_gradient_flow_dt_when_required(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_a = make_run(root, "run_a", gradient_dt=0.05)
            run_b = make_run(root, "run_b", seed=2, gradient_dt=0.05)
            run_c = make_run(root, "run_c", seed=3, gradient_dt=0.25)

            groups = finalize_analysis.group_run_directories(
                finalize_analysis.discover_run_directories(str(root)),
                require_fixed_dt=True,
            )

            grouped_paths = [paths for _, paths in groups]
            self.assertEqual(len(grouped_paths), 2)
            self.assertIn(sorted([str(run_a.resolve()), str(run_b.resolve())]), grouped_paths)
            self.assertIn([str(run_c.resolve())], grouped_paths)

    def test_run_evaluation_group_key_includes_gradient_flow_dt(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_a = make_run(root, "run_a", gradient_dt=0.05)
            run_b = make_run(root, "run_b", seed=2, gradient_dt=0.25)

            self.assertNotEqual(
                run_evaluation._group_key_for_run(str(run_a)),
                run_evaluation._group_key_for_run(str(run_b)),
            )

    def test_filter_run_directories_accepts_search_style_aliases(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_a = make_run(root, "run_a", beta=2.4, epsilon1=0.0)
            run_b = make_run(root, "run_b", beta=2.4, epsilon1=0.0, seed=2)
            make_run(root, "run_c", beta=2.6, epsilon1=0.0)

            criteria = finalize_analysis.parse_filter_tokens(["beta=2.4", "eps1=0.0", "L0=4"])
            matches = finalize_analysis.filter_run_directories(
                finalize_analysis.discover_run_directories(str(root)),
                criteria,
            )

            self.assertEqual(matches, sorted([str(run_a.resolve()), str(run_b.resolve())]))

    def test_exclude_run_directories_filters_path_text(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            keep = make_run(root / "20260515", "1")
            drop = make_run(root / "20260514_fast", "1", seed=2)

            filtered = finalize_analysis.exclude_run_directories(
                [str(keep), str(drop)],
                ["20260514_fast"],
            )

            self.assertEqual(filtered, [str(keep.resolve())])

    def test_custom_thermalization_threads_into_loader(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_a = make_run(root, "run_a")
            run_b = make_run(root, "run_b", seed=2, n_sweep=200, run_offset=0.005)

            combined_default, meta_default = run_evaluation._load_combined_w_temp(
                [str(run_a), str(run_b)],
                load_workers=1,
                thermalization_steps=1500,
            )
            combined_strict, meta_strict = run_evaluation._load_combined_w_temp(
                [str(run_a), str(run_b)],
                load_workers=1,
                thermalization_steps=2500,
            )

            self.assertIsNotNone(combined_default)
            self.assertIsNotNone(combined_strict)
            self.assertEqual(meta_default["n_configurations_after_cut"], 4)
            self.assertEqual(meta_strict["n_configurations_after_cut"], 2)

    def test_per_run_thermalization_threads_into_loader(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_a = make_run(root, "run_a")
            run_b = make_run(root, "run_b", seed=2, n_sweep=200, run_offset=0.005)

            combined, metadata = run_evaluation._load_combined_w_temp(
                [str(run_a), str(run_b)],
                load_workers=1,
                thermalization_steps_by_run={
                    str(run_a): 1500,
                    str(run_b): 2500,
                },
            )

            self.assertIsNotNone(combined)
            self.assertEqual(metadata["n_configurations_after_cut"], 3)
            self.assertIsNone(metadata["thermalization_steps"])
            self.assertEqual(
                metadata["thermalization_steps_by_run"],
                {
                    str(run_a.resolve()): 1500,
                    str(run_b.resolve()): 2500,
                },
            )

    def test_mixed_wilson_loop_merge_keeps_union_with_per_pair_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_a = make_run(root, "run_a", w_temp_r_values=[2, 3], w_temp_t_values=[1, 2])
            run_b = make_run(root, "run_b", seed=2, w_temp_r_values=[2, 3, 4], w_temp_t_values=[1, 2])

            combined, metadata = run_evaluation._load_combined_w_temp(
                [str(run_a), str(run_b)],
                load_workers=1,
                thermalization_steps=1500,
            )

            self.assertIsNotNone(combined)
            assert combined is not None
            self.assertIn((None, 4, 1), combined.wilson_by_flow_pair)
            self.assertEqual(len(combined.wilson_by_flow_pair[(None, 2, 1)]), 4)
            self.assertEqual(len(combined.wilson_by_flow_pair[(None, 4, 1)]), 2)
            self.assertEqual(metadata["wilson_loop_sample_counts"]["none,2,1"], 4)
            self.assertEqual(metadata["wilson_loop_sample_counts"]["none,4,1"], 2)
            calc = Calculator(combined, n_bootstrap=8, step_size=1)
            calc.prime_w_rt_cache(pairs=[(2, 1), (4, 1)])
            self.assertTrue(np.isfinite(calc.get_variable("W_R_T", R=4, T=1).get()))

    def test_gradient_flow_merge_keeps_union_with_per_flow_time_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_a = make_run(root, "run_a", gradient_t_values=[0.0, 1.0, 3.7])
            run_b = make_run(root, "run_b", seed=2, gradient_t_values=[0.0, 1.0, 3.56])
            write_gradient_flow_obs(
                run_a,
                [
                    (0, 0.0, 1.0, 0.00),
                    (0, 1.0, 2.0, 0.04),
                    (0, 3.7, 5.0, 0.12),
                    (1, 0.0, 1.1, 0.00),
                    (1, 1.0, 2.1, 0.05),
                    (1, 3.7, 5.1, 0.13),
                ],
            )
            write_gradient_flow_obs(
                run_b,
                [
                    (0, 0.0, 1.2, 0.00),
                    (0, 1.0, 2.2, 0.05),
                    (0, 3.56, 4.8, 0.11),
                    (1, 0.0, 1.3, 0.00),
                    (1, 1.0, 2.3, 0.06),
                    (1, 3.56, 4.9, 0.12),
                ],
            )

            combined, metadata = run_evaluation._load_combined_gradient_flow_obs(
                [str(run_a), str(run_b)],
                load_workers=1,
                thermalization_steps=0,
            )

            self.assertIsNotNone(combined)
            combined = combined or {}
            self.assertEqual(metadata["available_flow_times"], [0.0, 1.0, 3.56, 3.7])
            self.assertEqual(metadata["n_configurations_by_flow_time"]["1"], 4)
            self.assertEqual(metadata["n_configurations_by_flow_time"]["3.56"], 2)
            self.assertEqual(len(combined[3.56]["t2E_clover"]), 2)
            self.assertEqual(len(combined[3.7]["t2E_clover"]), 2)

    def test_gradient_flow_summary_keeps_unequal_flow_time_lengths(self):
        flow_data = {
            0.0: {
                "Ehat_clover": np.asarray([1.0, 2.0, 3.0], dtype=np.float32),
                "t2E_clover": np.asarray([0.1, 0.2, 0.3], dtype=np.float32),
            },
            0.5: {
                "Ehat_clover": np.asarray([2.0, 4.0], dtype=np.float32),
                "t2E_clover": np.asarray([0.2, 0.4], dtype=np.float32),
            },
        }

        summary = run_evaluation.summarize_gradient_flow_obs(
            flow_data,
            t0_target=0.3,
            block_size=1,
            n_bootstrap=8,
            include_bootstrap=True,
        )

        self.assertEqual(summary["n_configurations_used"], 2)
        self.assertEqual(summary["n_configurations_min"], 2)
        self.assertEqual(summary["n_configurations_max"], 3)
        self.assertFalse(summary["truncated_to_common_length"])
        self.assertAlmostEqual(summary["Ehat_clover"]["0"], 2.0)
        self.assertEqual(summary["bootstrap_samples"]["t2E_clover"].shape, (2, 8))

    def test_gradient_flow_summary_reports_fixed_0p1_crossing_and_weighted_fit(self):
        flow_data = {
            0.0: {
                "Ehat_clover": np.asarray([1.0, 1.0, 1.0], dtype=np.float32),
                "t2E_clover": np.asarray([0.040, 0.042, 0.044], dtype=np.float32),
            },
            0.5: {
                "Ehat_clover": np.asarray([1.0, 1.0, 1.0], dtype=np.float32),
                "t2E_clover": np.asarray([0.070, 0.072, 0.074], dtype=np.float32),
            },
            1.0: {
                "Ehat_clover": np.asarray([1.0, 1.0, 1.0], dtype=np.float32),
                "t2E_clover": np.asarray([0.120, 0.122, 0.124], dtype=np.float32),
            },
            1.5: {
                "Ehat_clover": np.asarray([1.0, 1.0, 1.0], dtype=np.float32),
                "t2E_clover": np.asarray([0.160, 0.162, 0.164], dtype=np.float32),
            },
        }

        summary = run_evaluation.summarize_gradient_flow_obs(
            flow_data,
            t0_target=0.3,
            block_size=1,
            n_bootstrap=8,
            include_bootstrap=True,
        )

        self.assertAlmostEqual(summary["t_over_a2_at_t2E_clover_0p1"], 0.78, places=6)
        self.assertIsNotNone(summary["t_over_a2_at_t2E_clover_0p1_weighted_fit"])
        self.assertEqual(len(summary["t_over_a2_at_t2E_clover_0p1_weighted_fit_points"]), 4)
        self.assertEqual(summary["bootstrap_samples"]["t_over_a2_at_t2E_clover_0p1"].shape, (8,))

    def test_whitespace_w_temp_loading(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "W_temp.out"
            path.write_text(
                "# step L T W_temp\n"
                "0 1 1 0.5\n"
                "0 1 2 0.25\n"
                "5 1 1 0.4\n"
                "5 1 2 0.2\n",
                encoding="utf-8",
            )

            fd = data_organizer.FileData(str(path)).read_file()
            self.assertEqual([obs.name for obs in fd.observables], ["step", "L", "T", "W_temp"])
            compact = data_organizer.load_compact_wilson_file(str(path), min_step=0)

            self.assertIsNotNone(compact)
            assert compact is not None
            self.assertEqual(compact.pair_order, [(1, 1), (1, 2)])
            self.assertEqual(compact.n_configurations, 2)

    def test_flowed_w_temp_loading_uses_flow_time_key(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "gradient_flow_wtemp.dat"
            path.write_text(
                "# conf_id t_over_a2 L T W_temp\n"
                "0 0.0 1 1 0.5\n"
                "0 0.25 1 1 0.6\n"
                "5 0.0 1 1 0.4\n"
                "5 0.25 1 1 0.55\n",
                encoding="utf-8",
            )

            compact = data_organizer.load_compact_wilson_file(str(path), min_step=0)

            self.assertIsNotNone(compact)
            assert compact is not None
            self.assertEqual(compact.available_flow_times(), [None] if False else [0.0, 0.25])
            calc = Calculator(compact, n_bootstrap=8, step_size=1)
            self.assertAlmostEqual(calc.get_variable("W_R_T", R=1, T=1, flow_time=0.25).get(), 0.575)

    def test_new_input_yaml_defaults_legacy_dimensions(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            payload = {
                "MetropolisParams": {
                    "L0": 4,
                    "L1": 4,
                    "L2": 4,
                    "L3": 4,
                    "nHits": 10,
                    "nSweep": 100,
                    "seed": 1,
                    "beta": 2.4,
                    "delta": 0.1,
                    "epsilon1": 0.0,
                    "epsilon2": 0.0,
                },
                "GaugeObservableParams": {"W_temp_L_T_pairs": [["1:2", "1:2"]]},
                "GradientFlowParams": {"enabled": True, "t_values": [0.0, 0.25]},
            }
            with (run_dir / "input.yaml").open("w", encoding="utf-8") as handle:
                yaml.safe_dump(payload, handle)

            metro, gauge = load_params(str(run_dir / "input.yaml"))

            self.assertEqual((metro.Ndims, metro.Nd, metro.Nc), (4, 4, 2))
            self.assertEqual(gauge.W_temp_L_T_pairs, [(1, 1), (1, 2), (2, 1), (2, 2)])

    def test_creutz_chi_uses_standard_adjacent_loop_definition(self):
        compact = data_organizer.CompactWilsonData(
            "synthetic",
            flow_pair_order=[(None, 1, 1), (None, 1, 2), (None, 2, 1), (None, 2, 2)],
            wilson_by_flow_pair={
                (None, 1, 1): np.asarray([0.5, 0.6], dtype=np.float32),
                (None, 1, 2): np.asarray([0.25, 0.35], dtype=np.float32),
                (None, 2, 1): np.asarray([0.2, 0.3], dtype=np.float32),
                (None, 2, 2): np.asarray([0.08, 0.12], dtype=np.float32),
            },
        )
        calc = Calculator(compact, n_bootstrap=8, step_size=1)

        chi = calc.get_variable("chi", R=1.5, T=1.5)
        loop_means = {
            pair: float(np.mean(values))
            for pair, values in compact.wilson_by_flow_pair.items()
        }
        expected = creutz_chi_from_wilson(
            loop_means[(None, 1, 1)],
            loop_means[(None, 1, 2)],
            loop_means[(None, 2, 1)],
            loop_means[(None, 2, 2)],
        )

        self.assertAlmostEqual(chi.get(), expected, places=7)

    def test_calculator_normalizes_equivalent_flow_time_cache_keys(self):
        compact = data_organizer.CompactWilsonData(
            "synthetic",
            flow_pair_order=[(None, 1, 1), (None, 1, 2), (None, 2, 1), (None, 2, 2)],
            wilson_by_flow_pair={
                (None, 1, 1): np.asarray([0.5, 0.6], dtype=np.float32),
                (None, 1, 2): np.asarray([0.25, 0.35], dtype=np.float32),
                (None, 2, 1): np.asarray([0.2, 0.3], dtype=np.float32),
                (None, 2, 2): np.asarray([0.08, 0.12], dtype=np.float32),
            },
        )
        calc = Calculator(compact, n_bootstrap=8, step_size=1)

        implicit = calc.get_variable("chi", R=1.5, T=1.5)
        explicit = calc.get_variable("chi", R=np.float64(1.5), T=1.5, flow_time=None)

        self.assertIs(implicit, explicit)

    def test_search_data_extracts_dynamic_creutz_fields(self):
        calc_data = {
            "chi": {"1.5,2.5": 0.12},
            "chi_err": {"1.5,2.5": 0.03},
            "F_chi": {"1.5": 0.22},
            "F_chi_err": {"1.5": 0.04},
        }

        self.assertEqual(search_data._extract_calculated_value("chi_R1p5_T2p5", calc_data), 0.12)
        self.assertEqual(search_data._extract_calculated_value("chi_R1p5_T2p5_err", calc_data), 0.03)
        self.assertEqual(search_data._extract_calculated_value("F_chi_R1p5", calc_data), 0.22)
        self.assertEqual(search_data._extract_calculated_value("F_chi_R1p5_err", calc_data), 0.04)

        with contextlib.redirect_stderr(io.StringIO()):
            criteria, outputs = search_data.parse_tokens(["chi_R1p5_T2p5=0.12", "F_chi_R1p5_err"])
        self.assertEqual(criteria, {"chi_R1p5_T2p5": 0.12})
        self.assertEqual(outputs, ["F_chi_R1p5_err"])
        self.assertTrue(
            search_data._matches_calculated_criteria(
                search_data.RowSpec("row", ".", None, False),
                criteria,
                calc_data,
            )
        )

    def test_finalized_search_summary_uses_only_finalized_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            analysis_dir = root / "beta_2.4__L_24x24x24x64__eps1_0__eps2_0__nrun_2__abc123"
            write_json(
                analysis_dir / "manifest.json",
                {
                    "analysis_id": "abc123",
                    "updated_at": "2026-05-17T00:00:00+00:00",
                    "status": {"derived_complete": True, "gradient_flow_complete": True, "creutz_complete": True},
                    "block_size": 400,
                    "thermalization_steps": 800,
                    "target_force": 1.65,
                },
            )
            write_json(
                analysis_dir / "input_runs.json",
                {
                    "run_dirs": [
                        "/home/pietl/masterarbeit/data/20260514/1",
                        "/home/pietl/masterarbeit/data/20260514/2",
                    ],
                    "compatibility_summary": {
                        "metropolis_common": {
                            "beta": 2.4,
                            "L0": 24,
                            "L1": 24,
                            "L2": 24,
                            "L3": 64,
                            "epsilon1": 0.0,
                            "epsilon2": 0.0,
                        }
                    },
                },
            )
            write_json(
                analysis_dir / "derived" / "summary.json",
                {
                    "r0": 4.0,
                    "r0_err": 0.1,
                    "a": 0.125,
                    "a_err": 0.003,
                    "length": 3.0,
                    "length_err": 0.2,
                    "epsilon_bar": 0.0,
                    "epsilon_bar_err": 0.0,
                },
            )
            write_json(
                analysis_dir / "gradient_flow" / "summary.json",
                {
                    "t_over_a2_at_t2E_clover_0p1": 0.96,
                    "t_over_a2_at_t2E_clover_0p1_err": 0.02,
                    "t_over_a2_at_t2E_clover_0p1_weighted_fit": 0.95,
                    "t_over_a2_at_t2E_clover_0p1_weighted_fit_err": 0.03,
                },
            )
            write_json(
                analysis_dir / "creutz" / "summary.json",
                {
                    "chi": {
                        "1.5,1.5": 0.20,
                        "2.5,2.5": 0.01,
                        "3.5,2.5": 0.50,
                    },
                    "chi_err": {
                        "1.5,1.5": 0.02,
                        "2.5,2.5": 0.02,
                        "3.5,2.5": 0.01,
                    },
                },
            )

            rows = search_data.search_data(root, mode="summary", criteria={"eps1": 0.0}, output=False)

            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertEqual(row["analysis_id"], "abc123")
            self.assertEqual(row["source_path"], "~/masterarbeit/data/20260514/1-2")
            self.assertEqual(row["L"], "24x24x24x64")
            self.assertEqual(row["n_runs"], 2)
            self.assertEqual(row["r0"], 4.0)
            self.assertEqual(row["r0_err"], 0.1)
            self.assertEqual(row["gf_t_over_a2"], 0.96)
            self.assertEqual(row["gf_t_over_a2_err"], 0.02)
            self.assertEqual(list(row["creutz_R_eq_T"]), ["1.5"])
            self.assertEqual(row["creutz_R_eq_T"]["1.5"]["err"], 0.02)

    def test_finalized_search_quick_filters_by_query(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = root / "beta_2.4__L_24x24x24x64__eps1_0__eps2_0__nrun_1__first"
            second = root / "beta_2.5__L_32x32x32x64__eps1_0__eps2_0__nrun_1__second"
            write_json(first / "manifest.json", {"analysis_id": "first", "status": {}})
            write_json(second / "manifest.json", {"analysis_id": "second", "status": {}})

            rows = search_data.search_data(root, mode="quick", query="second", output=False)

            self.assertEqual([row["analysis_id"] for row in rows], ["second"])

    def test_finalized_summary_can_omit_creutz_columns(self):
        rows = [{"analysis_id": "abc", "creutz_R_eq_T": {"1.5": {"value": 0.2, "err": 0.01}}}]

        columns = search_data._finalized_summary_columns(rows)
        self.assertIn("chi_R=T_1.5", columns)
        self.assertNotIn("chi_R=T_1.5", search_data.FINALIZED_SUMMARY_COLUMNS)
        self.assertIn("source_path", search_data.FINALIZED_SUMMARY_COLUMNS)
        self.assertNotIn("epsilon2", search_data.FINALIZED_SUMMARY_COLUMNS)
        self.assertNotIn("n_runs", search_data.FINALIZED_SUMMARY_COLUMNS)

    def test_finalized_table_rounds_values_to_uncertainties(self):
        row = {
            "r0": 4.3991963,
            "r0_err": 0.017930627,
            "a": 0.11360924,
            "a_err": 0.00046252487,
            "gf_t_over_a2": 0.968204,
            "gf_t_over_a2_err": 0.00014875497,
            "epsilon_bar": 0.0,
            "epsilon_bar_err": 0.0,
        }

        formatted = search_data._format_table_row(
            row,
            [
                "r0",
                "r0_err",
                "a",
                "a_err",
                "gf_t_over_a2",
                "gf_t_over_a2_err",
                "epsilon_bar",
                "epsilon_bar_err",
            ],
        )

        self.assertEqual(
            formatted,
            [
                "4.399",
                "0.018",
                "0.1136",
                "0.0005",
                "0.96820",
                "0.00015",
                "0",
                "0",
            ],
        )

    def test_l1_only_data_reports_no_standard_creutz_ratio(self):
        compact = data_organizer.CompactWilsonData(
            "synthetic",
            flow_pair_order=[(None, 1, 1), (None, 1, 2)],
            wilson_by_flow_pair={
                (None, 1, 1): np.asarray([0.5, 0.4], dtype=np.float32),
                (None, 1, 2): np.asarray([0.25, 0.2], dtype=np.float32),
            },
        )
        calc = Calculator(compact, n_bootstrap=8, step_size=1)

        with self.assertRaisesRegex(ValueError, "Not enough adjacent L values"):
            calc.get_variable("r0_chi")

    def test_r0_helper_matches_calculator_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_dir = make_run(root, "run_a")

            combined, _ = run_evaluation._load_combined_w_temp(
                [str(run_dir)],
                load_workers=1,
                thermalization_steps=1500,
            )
            self.assertIsNotNone(combined)
            calc = Calculator(combined, n_bootstrap=16, step_size=1)
            r0_var = calc.get_variable("r0", t_min=1, t_max=None, target_force=1.65, r_min=2)

            fit_rs = calc.get_unique_Rs(r_min=2)
            v_vars = [calc.get_variable("V_R", R=int(r), t_min=1, t_max=None) for r in fit_rs]
            fit_result = fit_r0_from_potential_data(
                np.asarray(fit_rs, dtype=float),
                np.asarray([var.get() for var in v_vars], dtype=float),
                errs=np.asarray([var.err() for var in v_vars], dtype=float),
                bootstrap_matrix=np.stack([var.bootstrap() for var in v_vars], axis=0),
                target_force=1.65,
                n_bootstrap=16,
            )

            self.assertAlmostEqual(r0_var.get(), fit_result["r0"], places=10)
            self.assertAlmostEqual(r0_var.err(), fit_result["r0_err"], places=10)

    def test_finalize_analysis_resumes_and_writes_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            data_root = root / "data"
            output_root = root / "finalized"
            run_a = make_run(data_root, "run_a")
            run_b = make_run(data_root, "run_b", seed=2, n_sweep=200, run_offset=0.005)

            patchers = [
                mock.patch.object(finalize_analysis, "save_thermalization_plot", side_effect=write_stub_plot),
                mock.patch.object(finalize_analysis, "save_bootstrap_block_size_plot", side_effect=write_stub_plot),
                mock.patch.object(finalize_analysis, "save_effective_mass_plot", side_effect=write_stub_plot),
                mock.patch.object(finalize_analysis, "save_r0_stability_plot", side_effect=write_stub_plot),
                mock.patch.object(finalize_analysis, "save_cornell_plot", side_effect=write_stub_plot),
                mock.patch.object(finalize_analysis, "save_gradient_flow_plot", side_effect=write_stub_plot),
                mock.patch.object(finalize_analysis, "save_creutz_plot", side_effect=write_stub_plot),
                mock.patch.object(finalize_analysis, "save_creutz_diagonal_plot", side_effect=write_stub_plot),
            ]
            for patcher in patchers:
                patcher.start()
            try:
                first_runner = finalize_analysis.FinalizedAnalysisRunner(
                    run_dirs=[str(run_a), str(run_b)],
                    output_root=output_root,
                    plot_mode="html",
                    load_workers=1,
                    calc_workers=None,
                    n_bootstrap=16,
                    target_force=1.65,
                    block_size_scan_values=[1, 2],
                    open_plots=False,
                    input_func=InputSequence(["", "", "1", "1"]),
                    print_func=lambda *args, **kwargs: None,
                )
                with self.assertRaisesRegex(RuntimeError, "input exhausted"):
                    first_runner.run()

                analysis_dir = first_runner.analysis_dir
                self.assertIsNotNone(first_runner.thermalization_preview_path)
                assert first_runner.thermalization_preview_path is not None
                self.assertTrue(first_runner.thermalization_preview_path.exists())
                first_choices = finalize_analysis.load_json(
                    analysis_dir / "wilsonloop" / "choices.json",
                    default={"choices": {}},
                )
                self.assertEqual(sorted(first_choices["choices"]), ["2"])
                block_choice = finalize_analysis.load_json(
                    analysis_dir / "bootstrap_block_size_choice.json",
                    default={},
                )
                self.assertEqual(block_choice["block_size"], 1)
                selection_state = finalize_analysis.load_json(
                    first_runner.selection_preview_dir / "thermalization_cuts.json",
                    default={},
                )
                self.assertEqual(
                    selection_state["cuts_by_run"],
                    {
                        str(run_a.resolve()): run_evaluation.THERMALIZATION_STEPS,
                        str(run_b.resolve()): run_evaluation.THERMALIZATION_STEPS,
                    },
                )

                second_messages = []
                second_runner = finalize_analysis.FinalizedAnalysisRunner(
                    run_dirs=[str(run_a), str(run_b)],
                    output_root=output_root,
                    plot_mode="html",
                    load_workers=1,
                    calc_workers=None,
                    n_bootstrap=16,
                    target_force=1.65,
                    block_size_scan_values=[1, 2],
                    open_plots=False,
                    input_func=InputSequence(["1", "1", "2"]),
                    print_func=lambda *args, **kwargs: second_messages.append(
                        " ".join(str(arg) for arg in args)
                    ),
                )
                second_runner.run()

                manifest = finalize_analysis.load_json(analysis_dir / "manifest.json")
                self.assertTrue(manifest["status"]["wilsonloop_complete"])
                self.assertTrue(manifest["status"]["r0_complete"])
                self.assertTrue(manifest["status"]["derived_complete"])
                self.assertTrue(manifest["status"]["gradient_flow_complete"])
                self.assertTrue(manifest["status"]["creutz_complete"])
                self.assertIn("result_summary_txt", manifest)
                self.assertIn("thermalization_preview_plot", manifest)
                self.assertEqual(manifest["block_size"], 1)
                self.assertEqual(second_runner.calc.step_size, 1)
                self.assertEqual(
                    manifest["thermalization_steps_by_run"],
                    {
                        str(run_a.resolve()): run_evaluation.THERMALIZATION_STEPS,
                        str(run_b.resolve()): run_evaluation.THERMALIZATION_STEPS,
                    },
                )

                self.assertTrue((analysis_dir / "scan_cache" / "wrt_scan.json").exists())
                self.assertTrue((analysis_dir / "scan_cache" / "bootstrap_block_size_scan.json").exists())
                self.assertTrue((analysis_dir / "wilsonloop" / "V_R_summary.json").exists())
                self.assertTrue((analysis_dir / "r0" / "r0_result.json").exists())
                self.assertTrue((analysis_dir / "derived" / "summary.json").exists())
                self.assertTrue((analysis_dir / "gradient_flow" / "summary.json").exists())
                self.assertTrue((analysis_dir / "creutz" / "summary.json").exists())
                summary_txt_path = analysis_dir / "result_summary.txt"
                self.assertTrue(summary_txt_path.exists())
                summary_txt = summary_txt_path.read_text(encoding="utf-8")
                self.assertIn("r0=", summary_txt)
                self.assertIn("length=", summary_txt)
                self.assertIn(" fm", summary_txt)
                self.assertIn("eps_bar=", summary_txt)
                self.assertIn("gf_t/a²=", summary_txt)
                self.assertTrue((analysis_dir / "derived" / "bootstrap" / "volume_r0.npy").exists())
                self.assertTrue((analysis_dir / "plots" / "bootstrap_block_size.html").exists())
                self.assertTrue((analysis_dir / "plots" / "thermalization_preview.html").exists())
                self.assertTrue((analysis_dir / "plots" / "creutz_ratios_R_eq_T.html").exists())
                self.assertTrue((analysis_dir / "plots" / "thermalization_by_run").exists())
                self.assertTrue(any("Result summary:" in message for message in second_messages))
                self.assertTrue(any("  r0 fit:" in message for message in second_messages))
                self.assertTrue(any("  Derived:" in message for message in second_messages))
            finally:
                for patcher in reversed(patchers):
                    patcher.stop()


if __name__ == "__main__":
    unittest.main()
