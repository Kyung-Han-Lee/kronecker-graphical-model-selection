import argparse
from replication.adni.figures.criterion_curves import draw_curves


def main():
    parser = argparse.ArgumentParser(description="Draw the ADNI criterion curves.")
    parser.add_argument("--paths", required=True)
    parser.add_argument("--verdicts", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    draw_curves(args.paths, args.verdicts, args.output)


if __name__ == "__main__":
    main()
