import argparse

from replication.adni.python.band_diagnostics import GROUPS, run


def main():
    parser = argparse.ArgumentParser(
        description="Compare ADNI ROI correlations across pass-band thirds."
    )
    parser.add_argument("--aggregated", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--groups", default=",".join(GROUPS))
    parser.add_argument("--splits", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20261101)
    args = parser.parse_args()
    run(
        args.aggregated,
        args.output,
        groups=args.groups.split(","),
        splits=args.splits,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
