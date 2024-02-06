import argparse
import csv
import re
from pathlib import Path
import typing as t

import anndata
import pandas as pd


def filter_crosswalk_table(table: pd.DataFrame, *columns: str) -> pd.DataFrame:
    """Filter the crosswalk table to only include specified columns.

    Also removes empty rows and cast values to string.

    Args:
        table (pd.DataFrame): Original full crosswalk table

    Returns:
        pd.DataFrame: Filtered table
    """
    return table[list(columns)].dropna().astype(str).drop_duplicates()


def generate_iri(label: str) -> str:
    """Create a temporary IRI based on a label.

    Args:
        label (str): Label for the row

    Returns:
        str: Temporary IRI
    """
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
    """Crosswalks the data adding CLIDs and match types using a crosswalk table.

    Args:
        matrix (anndata.AnnData): Data to crosswalk
        data_label_column (str): Column used to match against the table
        data_clid_column (str): Column to store CLIDs in
        data_match_column (str): Column to store match type in
        table (pd.DataFrame): Crosswalk table
        table_label_column (str): Column used to match against the data
        table_clid_column (str): Column storing CLIDs
        table_match_column (str): Column storing match type

    Returns:
        anndata.AnnData: Crosswalked data with CLIDs and match type added
    """
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

    _set_default_clid(merged_obs, data_clid_column, data_label_column)
    _set_default_match(merged_obs, data_match_column)

    result = matrix.copy()
    result.obs = merged_obs
    return result


def _set_default_clid(obs: pd.DataFrame, clid_column: str, label_column: str) -> None:
    """Adds default CLIDs to rows that did not match against the crosswalk table.

    Args:
        obs (pd.DataFrame): Data rows
        clid_column (str): Column to check and update with default CLIDs
        label_column (str): Column used when generating default CLIDs
    """
    defaults = obs.apply(lambda row: generate_iri(row[label_column]), axis=1)
    obs.loc[obs[clid_column].isna(), clid_column] = defaults


def _set_default_match(obs: pd.DataFrame, column: str) -> None:
    """Adds default match type to rows that did not match against the crosswalk table.

    Args:
        obs (pd.DataFrame): Data rows
        column (str): Column to check and update with default match type
    """
    obs.loc[obs[column].isna(), column] = "skos:exactMatch"


def _get_empty_table(args: argparse.Namespace) -> pd.DataFrame:
    """Creates an empty crosswalk table.

    Args:
        args (argparse.Namespace): Same arguments as provided to `main`

    Returns:
        pd.DataFrame: An empty table
    """
    return pd.DataFrame(
        columns=[
            args.crosswalk_table_label_column,
            args.crosswalk_table_clid_column,
            args.crosswalk_table_match_column,
        ]
    )


def _read_table(path: str) -> t.Optional[pd.DataFrame]:
    """Read a crosswalking table. Metadata rows before the header are skipped.

    Args:
        path (str): Path to the csv file

    Returns:
        pd.DataFrame: A data frame with the table data
    """
    with open(path) as file:
        for row in csv.reader(file):
            if row[0].lower() == 'organ_level':
                return pd.read_csv(file, names=row)
    return None


def main(args: argparse.Namespace):
    """Crosswalks a h5ad file and saves the result to another h5ad file.

    Args:
        args (argparse.Namespace):
            CLI arguments, must contain "matrix",
            "annotation_column", "clid_column", "match_column",
            "crosswalk_table", "crosswalk_table_label_column",
            "crosswalk_table_clid_column", "crosswalk_table_match_column", and
            "output_matrix"
    """
    matrix = crosswalk(
        args.matrix,
        args.annotation_column,
        args.clid_column,
        args.match_column,
        args.crosswalk_table
        if args.crosswalk_table is not None
        else _get_empty_table(args),
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
        type=_read_table,
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
