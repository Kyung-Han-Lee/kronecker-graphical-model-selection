from contextlib import contextmanager
import csv
import gzip
import json
import os
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd


def read_rows(path, required=()):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as stream:
        reader = csv.DictReader(stream)
        missing = set(required) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        return list(reader)


def read_frame(path, required=(), columns=None, precision="round_trip"):
    path = Path(path)
    header = pd.read_csv(path, nrows=0).columns
    missing = set(required) - set(header)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    dtype = {
        name: "category" for name in ("id", "criterion", "status") if name in header
    }
    return pd.read_csv(path, usecols=columns, dtype=dtype, float_precision=precision)


def finite(frame, columns):
    for name in columns:
        if not np.isfinite(pd.to_numeric(frame[name], errors="raise")).all():
            raise ValueError(f"nonfinite values in {name}")


def write_rows(path, rows):
    path = Path(path)
    rows = list(rows)
    if not rows:
        raise ValueError("no output rows")
    with path.open("w", newline="") as stream:
        keys = list(dict.fromkeys(key for row in rows for key in row))
        writer = csv.DictWriter(stream, fieldnames=keys, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, value):
    path = Path(path)
    if path.exists():
        raise FileExistsError(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(value, indent=2, allow_nan=False) + "\n"
    descriptor, temporary = tempfile.mkstemp(dir=path.parent, prefix=".partial-")
    try:
        with os.fdopen(descriptor, "w") as stream:
            stream.write(text)
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def check_output(path):
    path = Path(path)
    if path.exists():
        raise FileExistsError(path)
    if path.parent.exists() and not path.parent.is_dir():
        raise NotADirectoryError(path.parent)


def read_json(path):
    return json.loads(Path(path).read_text())


@contextmanager
def output_directory(path):
    path = Path(path)
    if path.exists():
        raise FileExistsError(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=path.parent, prefix=".partial-") as folder:
        staging = Path(folder) / "result"
        staging.mkdir()
        yield staging
        if path.exists():
            raise FileExistsError(path)
        staging.rename(path)
