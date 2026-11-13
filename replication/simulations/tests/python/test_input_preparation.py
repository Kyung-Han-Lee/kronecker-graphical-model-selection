import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[4]
SCRIPTS = ROOT / "replication/simulations/scripts"


class InputPreparationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.base = Path(cls.temporary.name)
        settings = cls.base / "design.csv"
        with settings.open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(
                (
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
            )
            writer.writerow(
                (
                    "toy_inputs",
                    2,
                    "er",
                    "3x3",
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
        for script, experiment in (
            ("run_matrix.R", "matrix"),
            ("run_refit_comparison.R", "refit_comparison"),
        ):
            cls.command(
                [
                    "Rscript",
                    str(SCRIPTS / script),
                    f"--settings={settings}",
                    f"--out={cls.base / 'raw' / experiment}",
                    "--save_inputs=1",
                    "--workers=1",
                ]
            )
        cls.inputs = cls.base / "inputs"
        cls.command(
            [
                "Rscript",
                str(SCRIPTS / "prepare_inputs.R"),
                f"--outputs={cls.base / 'raw'}",
                f"--out={cls.inputs}",
                "--profiles=all",
            ]
        )

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    @classmethod
    def command(cls, arguments):
        result = subprocess.run(arguments, cwd=ROOT, capture_output=True, text=True)
        if result.returncode:
            raise AssertionError(result.stdout + "\n" + result.stderr)
        return result

    def module(self, name, output, *arguments):
        return self.command(
            [
                sys.executable,
                "-m",
                f"replication.simulations.scripts.{name}",
                "--input",
                str(self.inputs),
                "--output",
                str(output),
                *map(str, arguments),
            ]
        )

    def test_exported_inputs_include_candidate_information(self):
        for name in (
            "paths.csv",
            "scores.csv",
            "supports.csv",
            "intervals.csv",
            "scales.csv",
            "conditions.csv",
            "covariances.csv",
            "cases.json",
        ):
            self.assertTrue((self.inputs / name).is_file(), name)
        with (self.inputs / "paths.csv").open() as stream:
            paths = list(csv.DictReader(stream))
        self.assertEqual(len(paths), 2)
        self.assertEqual(len({row["id"] for row in paths}), 2)
        cases = json.loads((self.inputs / "cases.json").read_text())
        self.assertEqual(len(cases), 2)
        for case in cases:
            payload = json.loads(
                (self.inputs / "profiles" / (case["id"] + ".json")).read_text()
            )
            self.assertEqual(len(payload["candidates"]), 40)

    def test_independent_scalar_commands(self):
        commands = (
            ("check_relaxation_intervals", ()),
            ("calculate_unit_multiplier", ()),
            ("check_scale_certificate", ("--dimensions", "3", "--sample-sizes", "20")),
            (
                "check_partial_covariances",
                ("--dimensions", "3", "--sample-sizes", "20"),
            ),
            ("check_condition_numbers", ("--dimensions", "3", "--sample-sizes", "20")),
        )
        for name, arguments in commands:
            with self.subTest(command=name):
                output = self.base / (name + "-result")
                self.module(name, output, *arguments)
                self.assertTrue(output.exists())

    def test_fixed_profile_command(self):
        output = self.base / "profile-result"
        self.module("check_profile_optimisation", output, "--max-boxes", "1")
        summary = json.loads((output / "summary.json").read_text())
        self.assertEqual(summary["paths"], 2)
        self.assertEqual(summary["candidates"], 80)

    def test_refit_comparison_and_convergence(self):
        output = self.base / "refit-result"
        self.module("compare_refit_likelihoods", output)
        jobs = output / "refit_jobs.csv"
        self.assertTrue(jobs.is_file())
        convergence = self.base / "convergence.csv"
        self.command(
            [
                "Rscript",
                str(SCRIPTS / "prepare_refit_convergence.R"),
                f"--outputs={self.base / 'raw'}",
                f"--inputs={self.inputs}",
                f"--jobs={jobs}",
                f"--out={convergence}",
            ]
        )
        result = self.base / "convergence-result.json"
        self.command(
            [
                sys.executable,
                "-m",
                "replication.simulations.scripts.check_refit_convergence",
                "--input",
                str(convergence),
                "--output",
                str(result),
            ]
        )
        self.assertGreater(json.loads(result.read_text())["candidates"], 0)

    def test_existing_output_is_rejected_before_reading_inputs(self):
        result = subprocess.run(
            [
                "Rscript",
                str(SCRIPTS / "prepare_inputs.R"),
                f"--outputs={self.base / 'missing'}",
                f"--out={self.inputs}",
            ],
            cwd=ROOT,
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("output already exists", result.stderr)


if __name__ == "__main__":
    unittest.main()
