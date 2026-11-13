import argparse
import json
from pathlib import Path
from ..python import refit


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(refit.run(args.input, args.output)))


if __name__ == "__main__":
    main()
