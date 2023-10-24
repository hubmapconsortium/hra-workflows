import argparse
import anndata
import pandas as pd
from pathlib import Path


def filter_crosswalk_table(
    crosswalk_table: pd.DataFrame,
    crosswalk_table_label_column: str,
    crosswalk_table_clid_column: str,
    crosswalk_table_match_column: str,
) -> pd.DataFrame:
    """Filters the table to remove empty rows and keep only necessary columns"""
    COLUMNS = [
        crosswalk_table_label_column,
        crosswalk_table_clid_column,
        crosswalk_table_match_column,
    ]
    crosswalk_table.dropna(inplace=True)
    crosswalk_table[COLUMNS] = crosswalk_table[COLUMNS].astype(str)
    return crosswalk_table[COLUMNS].drop_duplicates()


def crosswalk(
    matrix: anndata.AnnData,
    annotation_column: str,
    crosswalk_table: pd.DataFrame,
    crosswalk_table_label_column: str,
    crosswalk_table_clid_column: str,
    crosswalk_table_match_column: str,
) -> anndata.AnnData:
    """Gives each cell a CL ID and Match type using crosswalk table"""
    filtered_crosswalk_table = filter_crosswalk_table(
        crosswalk_table,
        crosswalk_table_label_column,
        crosswalk_table_clid_column,
        crosswalk_table_match_column,
    )
    merged_obs = matrix.obs.merge(
        filtered_crosswalk_table,
        left_on=annotation_column,
        right_on=crosswalk_table_label_column,
        how="left",
    ).drop(crosswalk_table_label_column, axis=1)
    merged_obs.index = matrix.obs.index
    matrix.obs = merged_obs
    return matrix


def main(args: argparse.Namespace):
    matrix = crosswalk(
        args.matrix,
        args.annotation_column,
        args.crosswalk_table,
        args.crosswalk_table_label_column,
        args.crosswalk_table_clid_column,
        args.crosswalk_table_match_column,
    )
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
        "--crosswalk-table-label-column",
        required=True,
        help="Column with Azimuth label in crosswalking table",
    )
    parser.add_argument(
        "--crosswalk-table-clid-column",
        required=True,
        help="Column with CL ID in crosswalking table",
    )
    parser.add_argument(
        "--crosswalk-table-match-column",
        required=True,
        help="Column with match type in crosswalking table",
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
