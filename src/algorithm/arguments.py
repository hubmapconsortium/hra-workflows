import argparse
import typing as t
from pathlib import Path


def add_common_arguments(
    parser: t.Optional[argparse.ArgumentParser] = None,
) -> argparse.ArgumentParser:
    """Add arguments common to all algorithms.

    Args:
        parser (argparse.ArgumentParser, optional): An existing parser to add argument to. Defaults to None.

    Returns:
        argparse.ArgumentParser: A parser with common arguments added
    """
    if parser is None:
        parser = argparse.ArgumentParser(description="Compute annotations")

    parser.add_argument("matrix", type=Path, help="h5ad data file")
    parser.add_argument("--organ", required=True, help="Organ uberon id")
    parser.add_argument(
        "--organ-metadata",
        type=Path,
        default="/organ-metadata.json",
        help="Organ metadata file",
    )
    parser.add_argument(
        "--output-matrix",
        type=Path,
        default="annotated_matrix.h5ad",
        help="Annotated matrix output path",
    )
    parser.add_argument(
        "--output-annotations",
        type=Path,
        default="annotations.csv",
        help="Annotations csv output path",
    )
    parser.add_argument(
        "--output-report",
        type=Path,
        default="report.json",
        help="Report json output path",
    )

    return parser
