import argparse
from replication.adni.figures.connectome import generate


def main():
    parser = argparse.ArgumentParser(description="Draw the ADNI ROI networks.")
    parser.add_argument("--fits", required=True)
    parser.add_argument("--atlas-image", required=True)
    parser.add_argument("--atlas-labels", required=True)
    parser.add_argument("--surface", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    generate(args.output, args.fits, args.atlas_image, args.atlas_labels, args.surface)


if __name__ == "__main__":
    main()
