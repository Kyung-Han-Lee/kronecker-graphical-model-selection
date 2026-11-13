import argparse
import json
from pathlib import Path
from ..python import profiles


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--fraction", type=float, default=0.01)
    parser.add_argument("--max-boxes", type=int, default=20000)
    args = parser.parse_args()
    print(
        json.dumps(
            profiles.run(
                args.input, args.output, args.workers, args.fraction, args.max_boxes
            )
        )
    )


if __name__ == "__main__":
    main()
