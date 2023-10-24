import argparse
import anndata
import pandas as pd
from pathlib import Path


def filter_crosswalk_table(crosswalk_table: pd.DataFrame) -> pd.DataFrame:
    """Filters the table to remove empty rows and keep only necessary columns"""
    crosswalk_table.dropna(inplace=True)
    crosswalk_table[["A_L", "CL_ID", "CL_Match"]] = crosswalk_table[
        ["A_L", "CL_ID", "CL_Match"]
    ].astype(str)
    return crosswalk_table[["A_L", "CL_ID", "CL_Match"]].drop_duplicates()


def crosswalk(
    matrix: anndata.AnnData, annotation_column: str, crosswalk_table: pd.DataFrame
) -> anndata.AnnData:
    """Gives each cell a CL ID and Match type using crosswalk table"""
    filtered_crosswalk_table = filter_crosswalk_table(crosswalk_table)
    merged_obs = matrix.obs.merge(
        filtered_crosswalk_table, left_on=annotation_column, right_on="A_L", how="left"
    ).drop("A_L", axis=1)
    merged_obs.index = matrix.obs.index
    matrix.obs = merged_obs
    return matrix


def main(args: argparse.Namespace):
    matrix = crosswalk(args.matrix, args.annotation_column, args.crosswalk_table)
    matrix.write_h5ad(args.output_matrix)


def _get_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Add crosswalking to h5ad data")
    parser.add_argument("matrix", type=anndata.read_h5ad, help="h5ad data file")
    parser.add_argument(
        "--annotation-column", required=True, help="Column with annotations"
    )
    parser.add_argument(
        "--crosswalk-table",
        type=pd.read_csv,
        required=True,
        default="all_labels.csv",
        help="crosswalking csv file path",
    )
    parser.add_argument(
        "--output-matrix",
        type=Path,
        default="matrix_with_crosswalking.h5ad",
        help="matrix with crosswalking output path",
    )

    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
