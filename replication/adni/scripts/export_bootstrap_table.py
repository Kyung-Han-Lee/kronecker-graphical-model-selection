import argparse
from replication.adni.python.bootstrap_table import run


def main():
    parser = argparse.ArgumentParser(description="Export the ADNI bootstrap table.")
    parser.add_argument("--summary", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    run(args.summary, args.metadata, args.output)


if __name__ == "__main__":
    main()
