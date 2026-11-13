import importlib
import tempfile
import unittest
from pathlib import Path


class PublicModuleTests(unittest.TestCase):
    def test_input_columns_are_checked_before_calculation(self):
        io = importlib.import_module("replication.io")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "scores.csv"
            path.write_text("id,k\na,1\n")
            with self.assertRaisesRegex(ValueError, "q"):
                io.read_rows(path, required=("id", "k", "q"))

    def test_existing_output_directory_is_rejected_without_writes(self):
        io = importlib.import_module("replication.io")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "existing"
            path.mkdir()
            (path / "keep.txt").write_text("unchanged")
            with self.assertRaises(FileExistsError):
                with io.output_directory(path):
                    self.fail("existing output must be rejected before calculation")
            self.assertEqual((path / "keep.txt").read_text(), "unchanged")

    def test_failed_output_does_not_leave_partial_results(self):
        io = importlib.import_module("replication.io")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "result"
            with self.assertRaisesRegex(ValueError, "invalid input"):
                with io.output_directory(path) as temporary:
                    (temporary / "partial.csv").write_text("bad")
                    raise ValueError("invalid input")
            self.assertFalse(path.exists())

    def test_interval_calls_do_not_share_input_state(self):
        intervals = importlib.import_module("replication.simulations.python.intervals")
        paths = [{"id": "a", "experiment": "matrix", "M": "2"}]
        first = [{"id": "a", "k": "1", "q": "0"}, {"id": "a", "k": "2", "q": "2"}]
        second = [{"id": "a", "k": "1", "q": "0"}, {"id": "a", "k": "2", "q": "1"}]
        self.assertEqual(intervals.summarize(paths, first)["violations"], 0)
        self.assertEqual(intervals.summarize(paths, second)["violations"], 1)
        self.assertEqual(intervals.summarize(paths, first)["violations"], 0)


if __name__ == "__main__":
    unittest.main()
