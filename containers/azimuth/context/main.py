import argparse
import itertools
import json
import os
import subprocess
from pathlib import Path

import scanpy

# From azimuth-annotate docker image
import write_metadata


# Preferably theses values should be read from somewhere but they are
# currently only embedded in the azimuth_analysis R script
REFERENCES = ["RK", "LK", "RL", "LL", "HT"]


def _organ_or_reference(value: str):
    with open("/organ-mapping.json") as mapping_file:
        mapping = json.load(mapping_file)

    value = value.lower()
    items = itertools.chain(mapping.items(), zip(REFERENCES, REFERENCES))
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
        "-o", "--output", type=Path, default="annotations.csv", help="Output file"
    )
    parser.add_argument(
        "--output-matrix",
        type=Path,
        default="annotated_matrix.h5ad",
        help="Annotated matrix output file",
    )

    return parser


def main(args: argparse.Namespace):
    annotate_outputs = (
        "./secondary_analysis.h5ad",
        "./version_metadata.json",
        "./annotations.csv",
    )
    subprocess.run(
        [
            "Rscript",
            "/azimuth_analysis.R",
            args.reference,
            args.data,
            args.data,
        ],
        capture_output=True,
        check=True,
    )

    write_metadata.main(*annotate_outputs)
    scanpy.read_h5ad(annotate_outputs[0]).obs.to_csv(args.output)
    os.rename(annotate_outputs[0], args.output_matrix)


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
