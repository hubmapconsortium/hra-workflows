import argparse
import itertools
import json
import subprocess
from pathlib import Path

# Preferably theses values should be read from somewhere but they are
# currently only embedded in the azimuth_analysis R script
REFERENCES = ["RK", "LK", "RL", "LL", "HT"]


def _organ_or_reference(value: str):
    with open("./organ-mapping.json") as mapping_file:
        mapping = json.load(mapping_file)

    value = value.lower()
    items = itertools.chain(mapping.values(), zip(REFERENCES, REFERENCES))
    for key, reference in items:
        if key.lower() == value:
            return reference

    raise ValueError("Invalid organ")


def _get_arg_parser():
    parser = argparse.ArgumentParser(description="Compute annotations using azimuth")
    parser.add_argument("data", type=Path, help="h5ad data file")
    parser.add_argument(
        "--organ",
        type=_organ_or_reference,
        required=True,
        dest="reference",
        help="Organ uberon id",
    )
    parser.add_argument(
        "-o", "--output", type=Path, help="Output file", default="annotations.csv"
    )

    return parser


def main(args: argparse.Namespace):
    command = [
        "Rscript",
        "/azimuth_analysis.R",
        args.reference,
        args.data,
        args.data,
        args.output,
    ]
    subprocess.run(command, capture_output=True, check=True)


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
