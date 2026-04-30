from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import yaml

import finalize_analysis
import finalized_analysis_helpers as helpers
import run_evaluation
from calculator import Calculator, fit_r0_from_potential_data


def write_input_yaml(
    run_dir: Path,
    *,
    beta: float = 2.4,
    epsilon1: float = 0.2,
    epsilon2: float = 0.0,
    seed: int = 1,
    n_sweep: int = 100,
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
            "W_temp_L_T_pairs": [[2, "1:4"], [3, "1:4"], [4, "1:4"]],
            "W_mu_nu_pairs": [],
            "W_Lmu_Lnu_pairs": [],
            "plaquette_filename": "plaquette.out",
            "W_temp_filename": "W_temp.out",
            "W_mu_nu_filename": "W_mu_nu.out",
            "RetraceU_filename": "RetraceU.out",
            "write_to_file": True,
        },
    }
    with (run_dir / "input.yaml").open("w", encoding="utf-8") as handle:
        yaml.safe_dump(payload, handle, sort_keys=False)


def write_w_temp(run_dir: Path, *, run_offset: float = 0.0) -> None:
    A = 0.5
    sigma = 0.2
    B = 0.1
    steps = [0, 1000, 2000, 3000]
    rows = ["step,L,T,W_temp"]
    for cfg_index, step in enumerate(steps):
        for r_val in [2, 3, 4]:
            potential = A + sigma * r_val - B / r_val
            for t_val in [1, 2, 3, 4]:
                modulation = 1.0 + 0.01 * cfg_index + run_offset
                value = modulation * np.exp(-potential * t_val)
                rows.append(f"{step},{r_val},{t_val},{value:.10f}")
    (run_dir / "W_temp.out").write_text("\n".join(rows) + "\n", encoding="utf-8")


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
    )
    write_w_temp(run_dir, run_offset=run_offset)
    return run_dir


def write_stub_plot(path: Path, *args, **kwargs) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("<html><body>stub</body></html>", encoding="utf-8")


class InputSequence:
    def __init__(self, responses):
        self._responses = iter(responses)

    def __call__(self, prompt: str) -> str:
        try:
            return next(self._responses)
        except StopIteration as exc:
            raise RuntimeError("input exhausted") from exc


class FinalizeAnalysisTests(unittest.TestCase):
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
                    print_func=lambda *args, **kwargs: None,
                )
                second_runner.run()

                manifest = finalize_analysis.load_json(analysis_dir / "manifest.json")
                self.assertTrue(manifest["status"]["wilsonloop_complete"])
                self.assertTrue(manifest["status"]["r0_complete"])
                self.assertTrue(manifest["status"]["derived_complete"])
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
                self.assertTrue((analysis_dir / "derived" / "bootstrap" / "volume_r0.npy").exists())
                self.assertTrue((analysis_dir / "plots" / "bootstrap_block_size.html").exists())
                self.assertTrue((analysis_dir / "plots" / "thermalization_preview.html").exists())
                self.assertTrue((analysis_dir / "plots" / "thermalization_by_run").exists())
            finally:
                for patcher in reversed(patchers):
                    patcher.stop()


if __name__ == "__main__":
    unittest.main()
