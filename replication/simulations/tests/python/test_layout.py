from pathlib import Path
import csv
import shutil
import subprocess
import sys
import tempfile
import unittest

SIMULATIONS = Path(__file__).resolve().parents[2]


class SimulationLayoutTest(unittest.TestCase):
    def test_independent_commands_show_help_from_the_repository_root(self):
        names = (
            "check_relaxation_intervals.py",
            "calculate_unit_multiplier.py",
            "check_scale_certificate.py",
            "check_partial_covariances.py",
            "check_condition_numbers.py",
            "check_profile_optimisation.py",
            "compare_refit_likelihoods.py",
            "check_refit_convergence.py",
        )
        with tempfile.TemporaryDirectory():
            for name in names:
                with self.subTest(script=name):
                    result = subprocess.run(
                        [
                            sys.executable,
                            "-m",
                            "replication.simulations.scripts." + Path(name).stem,
                            "--help",
                        ],
                        cwd=SIMULATIONS.parents[1],
                        capture_output=True,
                        text=True,
                    )
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertIn("usage:", result.stdout)

    def test_matrix_command_writes_a_complete_small_dataset(self):
        fields = (
            "setting",
            "M",
            "family",
            "dims",
            "d",
            "m_new",
            "n",
            "n_idx",
            "reps",
            "w",
            "w_spread",
            "margin",
            "gradient",
        )
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            settings = base / "matrix.csv"
            with settings.open("w", newline="") as stream:
                writer = csv.writer(stream)
                writer.writerow(fields)
                writer.writerow(
                    (
                        "cli_matrix",
                        2,
                        "er",
                        "3x2",
                        1,
                        "",
                        20,
                        1,
                        1,
                        0.3,
                        0.25,
                        0.5,
                        "numeric",
                    )
                )
            output = base / "matrix"
            result = subprocess.run(
                [
                    "Rscript",
                    str(SIMULATIONS / "scripts/run_matrix.R"),
                    f"--settings={settings}",
                    f"--out={output}",
                    "--workers=1",
                ],
                cwd=directory,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("1 done, 0 reused, 0 failed", result.stdout)
            self.assertTrue((output / "datasets/cli_matrix_r001.rds").is_file())

    def test_stored_results_regenerate_tables_and_figures(self):
        stored = SIMULATIONS / "results"
        with tempfile.TemporaryDirectory() as directory:
            results = Path(directory) / "review"
            shutil.copytree(stored / "per_dataset", results / "per_dataset")
            subprocess.run(
                [
                    "Rscript",
                    str(SIMULATIONS / "scripts/summarize_results.R"),
                    "--from=per_dataset",
                    f"--results={results}",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            for folder, pattern in (("tables", "*.csv"), ("figures", "*.csv")):
                expected = sorted((stored / folder).glob(pattern))
                self.assertEqual(len(expected), 8 if folder == "tables" else 4)
                for path in expected:
                    self.assertEqual(
                        path.read_bytes(),
                        (results / folder / path.name).read_bytes(),
                        str(path),
                    )
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "replication.simulations.figures.draw_figures",
                    f"--results={results}",
                ],
                check=True,
                capture_output=True,
                text=True,
                cwd=SIMULATIONS.parents[1],
            )
            self.assertEqual(
                sorted(p.name for p in (results / "figures").glob("*.pdf")),
                ["Figure1.pdf", "Figure2.pdf", "Figure3.pdf"],
            )


if __name__ == "__main__":
    unittest.main()
