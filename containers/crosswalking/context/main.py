import argparse
import csv
import re
from pathlib import Path
import typing as t

import anndata
import pandas as pd


def filter_crosswalk_table(
    table: pd.DataFrame,
    organ_id: str,
    organ_level: str,
    organ_id_column: str,
    organ_level_column: str,
    table_label_column: str,
) -> pd.DataFrame:
    """Filter the crosswalk table to only include rows with organ id and level.

    Also removes empty rows and cast values to string.

    Args:
        table (pd.DataFrame): Original full crosswalk table

    Returns:
        pd.DataFrame: Filtered table
    """
    organ_id_rows = table[organ_id_column].str.lower() == organ_id.lower()
    organ_level_rows = table[organ_level_column].str.lower() == organ_level.lower()
    filtered_table = table[organ_id_rows & organ_level_rows]
    normalized_table = filtered_table.dropna().astype(str)
    unique_table = normalized_table.drop_duplicates(table_label_column)
    return unique_table


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
    organ_id: str,
    organ_level: str,
    data_label_column: str,
    data_clid_column: str,
    data_match_column: str,
    table: pd.DataFrame,
    table_organ_id_column: str,
    table_organ_level_column: str,
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
        table,
        organ_id,
        organ_level,
        table_organ_id_column,
        table_organ_level_column,
        table_label_column,
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
            args.crosswalk_table_organ_id_column,
            args.crosswalk_table_organ_level_column,
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
            if row[0].lower() == "organ_level":
                return pd.read_csv(file, names=row)
    return None


def main(args: argparse.Namespace):
    """Crosswalks a h5ad file and saves the result to another h5ad file.

    Args:
        args (argparse.Namespace):
            CLI arguments, must contain "matrix",
            "annotation_column", "clid_column", "match_column",
            "crosswalk_table", "crosswalk_table_organ_id_column",
            "crosswalk_table_organ_level_column", "crosswalk_table_label_column",
            "crosswalk_table_clid_column", "crosswalk_table_match_column", and
            "output_matrix"
    """
    metadata = args.matrix.uns["hra_crosswalking"]
    matrix = crosswalk(
        args.matrix,
        metadata["organ_id"],
        metadata["organ_level"],
        args.annotation_column,
        args.clid_column,
        args.match_column,
        (
            args.crosswalk_table
            if args.crosswalk_table is not None
            else _get_empty_table(args)
        ),
        args.crosswalk_table_organ_id_column,
        args.crosswalk_table_organ_level_column,
        args.crosswalk_table_label_column,
        args.crosswalk_table_clid_column,
        args.crosswalk_table_match_column,
    )
    matrix.write_h5ad(args.output_matrix)


def _get_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Add crosswalking to h5ad data")
    parser.add_argument("matrix", type=anndata.read_h5ad, help="h5ad data file")
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
        "--crosswalk-table",
        type=_read_table,
        help="crosswalking csv file path",
    )
    parser.add_argument(
        "--crosswalk-table-organ-id-column",
        default="Organ_ID",
        help="Column with organ ids in crosswalking table",
    )
    parser.add_argument(
        "--crosswalk-table-organ-level-column",
        default="Organ_Level",
        help="Column with organ levels in crosswalking table",
    )
    parser.add_argument(
        "--crosswalk-table-label-column",
        default="Annotation_Label",
        help="Column with label in crosswalking table",
    )
    parser.add_argument(
        "--crosswalk-table-clid-column",
        default="CL_ID",
        help="Column with CL ID in crosswalking table",
    )
    parser.add_argument(
        "--crosswalk-table-match-column",
        default="CL_Match",
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
