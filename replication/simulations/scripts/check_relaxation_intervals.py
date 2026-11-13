import argparse
import json
from pathlib import Path
from replication import io
from ..python import intervals
from ..python.inputs import load_paths


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    io.check_output(args.output)
    result = intervals.summarize(
        load_paths(args.input),
        io.read_frame(args.input / "intervals.csv", required=("id", "k", "q")),
    )
    io.write_json(args.output, result)
    print(json.dumps(result))


if __name__ == "__main__":
    main()
