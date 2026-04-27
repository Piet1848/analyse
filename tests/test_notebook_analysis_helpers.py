from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np
import yaml

import notebook_analysis_helpers as nah
import run_evaluation
from calculator import Calculator, cornell_potential_ansatz


def write_input_yaml(run_dir: Path, *, seed: int = 1, n_sweep: int = 100) -> None:
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
            "beta": 2.4,
            "delta": 0.1,
            "epsilon1": 0.2,
            "epsilon2": 0.0,
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
    a_param = 0.5
    sigma = 0.2
    b_param = 0.1
    steps = [0, 1000, 2000, 3000, 4000]
    rows = ["step,L,T,W_temp"]
    for cfg_index, step in enumerate(steps):
        for r_val in [2, 3, 4]:
            potential = a_param + sigma * r_val - b_param / r_val
            for t_val in [1, 2, 3, 4]:
                modulation = 1.0 + 0.01 * cfg_index + run_offset
                value = modulation * np.exp(-potential * t_val)
                rows.append(f"{step},{r_val},{t_val},{value:.10f}")
    (run_dir / "W_temp.out").write_text("\n".join(rows) + "\n", encoding="utf-8")


def make_run(base_dir: Path, name: str, *, seed: int = 1, n_sweep: int = 100, run_offset: float = 0.0) -> Path:
    run_dir = base_dir / name
    run_dir.mkdir(parents=True, exist_ok=True)
    write_input_yaml(run_dir, seed=seed, n_sweep=n_sweep)
    write_w_temp(run_dir, run_offset=run_offset)
    return run_dir


class NotebookAnalysisHelperTests(unittest.TestCase):
    def test_cache_key_stability_and_staleness_explanation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_a = make_run(root, "run_a")
            cuts = nah.resolve_thermalization_steps_by_run([str(run_a)], thermalization_steps=1500)
            options = nah.build_analysis_options(thermalization_steps_by_run=cuts, n_bootstrap=8)

            input_hash = nah.input_hash_for_runs([str(run_a)])
            key_a = nah.build_scan_cache_key([str(run_a)], options, input_hash=input_hash)
            key_b = nah.build_scan_cache_key([str(run_a)], options, input_hash=input_hash)

            self.assertEqual(key_a, key_b)

            expected = nah.expected_cache_metadata([str(run_a)], options, input_hash=input_hash)
            valid_payload = {
                **expected,
                "wrt_scan": [],
                "effective_mass_scan": [],
                "v_r_scan": [],
            }
            self.assertEqual(nah.explain_scan_cache_staleness(valid_payload, expected), [])

            stale_payload = dict(valid_payload)
            stale_payload["calculator_version"] = "old"
            self.assertIn("calculator_version differs", nah.explain_scan_cache_staleness(stale_payload, expected))

    def test_primed_wrt_cache_matches_lazy_calculation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_a = make_run(root, "run_a")
            combined, _ = run_evaluation._load_combined_w_temp(
                [str(run_a)],
                load_workers=1,
                thermalization_steps=1500,
            )
            self.assertIsNotNone(combined)
            assert combined is not None

            lazy = Calculator(combined, n_bootstrap=16, step_size=1)
            primed = Calculator(combined, n_bootstrap=16, step_size=1)
            primed.prime_w_rt_cache([(2, 1)])

            lazy_var = lazy.get_variable("W_R_T", R=2, T=1)
            primed_var = primed.get_variable("W_R_T", R=2, T=1)

            self.assertAlmostEqual(lazy_var.get(), primed_var.get(), places=7)
            self.assertAlmostEqual(lazy_var.err(), primed_var.err(), places=7)
            np.testing.assert_allclose(lazy_var.bootstrap(), primed_var.bootstrap(), rtol=1e-6)

    def test_physics_audit_reports_synthetic_cornell_and_creutz_checks(self):
        target_force = 1.65
        sigma = 0.2
        b_param = 0.1
        r0_val = float(np.sqrt((target_force - b_param) / sigma))
        wrt_records = []
        for r_val in [1, 2, 3]:
            for t_val in [1, 2, 3]:
                value = float(np.exp(-sigma * r_val * t_val))
                wrt_records.append({"R": r_val, "T": t_val, "value": value, "err": value * 0.01})

        audit = nah.build_physics_audit(
            {
                "analysis_options": {"target_force": target_force},
                "wrt_scan": wrt_records,
                "v_r_scan": [],
                "r0_scan": [
                    {
                        "r_min": 1,
                        "r0": r0_val,
                        "err": 0.1,
                        "A": 0.5,
                        "sigma": sigma,
                        "B": b_param,
                        "chi2_dof": 1.0,
                    }
                ],
                "bootstrap_quality": [
                    {
                        "variable": "W_R_T",
                        "params": {"R": 1, "T": 1},
                        "n_bootstrap": 8,
                        "finite_fraction": 1.0,
                        "invalid_count_before_repair": 0,
                        "error": 0.01,
                    }
                ],
            }
        )

        cornell_check = audit["cornell_sommer"]["checks"][0]
        self.assertAlmostEqual(cornell_check["target_force_delta"], 0.0, places=12)
        self.assertGreater(audit["creutz_ratio"]["checks"][0]["median_current_chi"], 0.0)

    def test_audit_does_not_change_current_formula_outputs(self):
        before = cornell_potential_ansatz(2.0, 0.5, 0.2, 0.1)
        audit = nah.build_physics_audit(
            {
                "analysis_options": {"target_force": 1.65},
                "wrt_scan": [],
                "v_r_scan": [],
                "r0_scan": [
                    {
                        "r_min": 1,
                        "r0": 2.0,
                        "err": 0.1,
                        "A": 0.5,
                        "sigma": -0.2,
                        "B": 2.0,
                        "chi2_dof": None,
                    }
                ],
                "bootstrap_quality": [],
            }
        )
        after = cornell_potential_ansatz(2.0, 0.5, 0.2, 0.1)

        self.assertEqual(before, after)
        self.assertGreater(audit["summary"]["warning_count"], 0)

    def test_load_or_build_scan_cache_writes_metadata_and_reuses_cache(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cache_root = root / "cache"
            run_a = make_run(root, "run_a")
            run_b = make_run(root, "run_b", seed=2, n_sweep=200, run_offset=0.005)

            first = nah.load_or_build_scan_cache(
                [str(run_a), str(run_b)],
                cache_root=cache_root,
                thermalization_steps=1500,
                n_bootstrap=8,
                block_size=1,
                load_workers=1,
            )
            second = nah.load_or_build_scan_cache(
                [str(run_a), str(run_b)],
                cache_root=cache_root,
                thermalization_steps=1500,
                n_bootstrap=8,
                block_size=1,
                load_workers=1,
            )

            self.assertFalse(first["_cache_status"]["loaded_from_cache"])
            self.assertTrue(second["_cache_status"]["loaded_from_cache"])
            self.assertEqual(first["cache_schema_version"], nah.CACHE_SCHEMA_VERSION)
            self.assertEqual(first["source_run_dirs"], [str(run_a.resolve()), str(run_b.resolve())])
            self.assertIn("audit", first)
            self.assertIn("wrt_scan", first)
            self.assertIn("effective_mass_scan", first)
            self.assertIn("v_r_scan", first)


if __name__ == "__main__":
    unittest.main()
