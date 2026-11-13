import argparse
import json
from pathlib import Path
from ..python import convergence
from replication import io


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    io.check_output(args.output)
    frame = io.read_frame(args.input, required=("saved_refit_ll", "strict_refit_ll"))
    result = convergence.summarize(frame)
    io.write_json(args.output, result)
    print(json.dumps(result))


if __name__ == "__main__":
    main()
