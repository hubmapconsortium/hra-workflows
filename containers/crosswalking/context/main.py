import argparse
import re
from pathlib import Path

import anndata
import pandas as pd


def filter_crosswalk_table(table: pd.DataFrame, *columns: str) -> pd.DataFrame:
    """Filters the table to remove empty rows and keep only necessary columns"""
    return table[list(columns)].dropna().astype(str).drop_duplicates()


def generate_iri(label: str):
    """generate IRIs for labels not found in crosswalk tables"""
    suffix = label.lower().strip()
    suffix = re.sub(r"\W+", "-", suffix)
    suffix = re.sub(r"[^a-z0-9-]+", "", suffix)
    return "ASCTB-TEMP:" + suffix


def crosswalk(
    matrix: anndata.AnnData,
    data_label_column: str,
    data_clid_column: str,
    data_match_column: str,
    table: pd.DataFrame,
    table_label_column: str,
    table_clid_column: str,
    table_match_column: str,
) -> anndata.AnnData:
    """Gives each cell a CL ID and Match type using crosswalk table"""
    column_map = {
        table_clid_column: data_clid_column,
        table_match_column: data_match_column,
    }
    table = filter_crosswalk_table(
        table, table_label_column, table_clid_column, table_match_column
    )
    merged_obs = (
        matrix.obs.merge(
            table, left_on=data_label_column, right_on=table_label_column, how="left"
        )
        .drop(columns=table_label_column)
        .rename(columns=column_map)
    )
    merged_obs.index = matrix.obs.index

    _set_default_clid(merged_obs, table_clid_column, data_label_column)
    _set_default_match(merged_obs, data_match_column)

    result = matrix.copy()
    result.obs = merged_obs
    return result


def _set_default_clid(obs: pd.DataFrame, clid_column: str, label_column: str):
    defaults = obs.apply(lambda row: generate_iri(row[label_column]), axis=1)
    obs.loc[obs[clid_column].isna(), clid_column] = defaults


def _set_default_match(obs: pd.DataFrame, column: str):
    obs.loc[obs[column].isna(), column] = "skos:exactMatch"


def _get_empty_table(args: argparse.Namespace) -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            args.crosswalk_table_label_column,
            args.crosswalk_table_clid_column,
            args.crosswalk_table_match_column,
        ]
    )


def main(args: argparse.Namespace):
    matrix = crosswalk(
        args.matrix,
        args.annotation_column,
        args.clid_column,
        args.match_column,
        args.crosswalk_table if args.crosswalk_table is not None else _get_empty_table(args),
        args.crosswalk_table_label_column,
        args.crosswalk_table_clid_column,
        args.crosswalk_table_match_column,
    )
    matrix.write_h5ad(args.output_matrix)


def _get_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Add crosswalking to h5ad data")
    parser.add_argument("matrix", type=anndata.read_h5ad, help="h5ad data file")
    parser.add_argument(
        "--crosswalk-table",
        type=pd.read_csv,
        help="crosswalking csv file path",
    )
    parser.add_argument(
        "--crosswalk-table-label-column",
        default="label",
        help="Column with Azimuth label in crosswalking table",
    )
    parser.add_argument(
        "--crosswalk-table-clid-column",
        default="clid",
        help="Column with CL ID in crosswalking table",
    )
    parser.add_argument(
        "--crosswalk-table-match-column",
        default="match",
        help="Column with match type in crosswalking table",
    )
    parser.add_argument(
        "--annotation-column", default="hra_prediction", help="Column with annotations"
    )
    parser.add_argument(
        "--clid-column", default="clid", help="Output column for cell ids"
    )
    parser.add_argument(
        "--match-column", default="match_type", help="Output column for match"
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
