import argparse
from replication.adni.python import projection


def main():
    parser = argparse.ArgumentParser(
        description="Project preprocessed time series onto the passband basis."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--groups", default="CN,MCI,Dementia")
    parser.add_argument("--option", default="C")
    args = parser.parse_args()
    projection.run(args.input, args.output, args.groups.split(","), args.option)


if __name__ == "__main__":
    main()
