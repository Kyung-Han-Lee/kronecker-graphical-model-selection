import argparse
import json
from pathlib import Path
from replication import io
from ..python import covariances
from ..python.inputs import load_paths


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dimensions", default="40,60,80")
    parser.add_argument("--sample-sizes", default="100,200,800,3200")
    args = parser.parse_args()
    io.check_output(args.output)
    dimensions = tuple(map(int, args.dimensions.split(",")))
    sample_sizes = tuple(map(int, args.sample_sizes.split(",")))
    paths = load_paths(args.input)
    values = io.read_frame(
        args.input / "covariances.csv", required=["id", "k", "minimum_eigenvalue"]
    )
    result = covariances.summarize(paths, values, dimensions, sample_sizes)
    io.write_json(args.output, result)
    print(json.dumps(result))


if __name__ == "__main__":
    main()
