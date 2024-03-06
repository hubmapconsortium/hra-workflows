import argparse
import csv
import typing as t
from pathlib import Path

import anndata
import pandas as pd

from src.util.ids import create_temp_asctb_id


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
    table_clid_label_column: str,
    table_match_column: str,
) -> anndata.AnnData:
    """Crosswalks the data adding CLIDs and match types using a crosswalk table.

    Args:
        matrix (anndata.AnnData): Data to crosswalk
        organ_id (str): Organ id
        organ_label (str): Organ level
        data_label_column (str): Column used to match against the table
        data_clid_column (str): Column to store CLIDs in
        data_match_column (str): Column to store match type in
        table (pd.DataFrame): Crosswalk table
        table_organ_id_column (str): Column storing organ uberon ids
        table_organ_level_column (str): Column storing organ levels
        table_label_column (str): Column used to match against the data
        table_clid_column (str): Column storing CLIDs
        table_clid_label_column (str): Column storing CL labels
        table_match_column (str): Column storing match type

    Returns:
        anndata.AnnData: Crosswalked data with CLIDs and match type added
    """
    column_map = {
        table_clid_column: data_clid_column,
        table_match_column: data_match_column,
    }
    matrix = _filter_invalid_rows(matrix, data_label_column)
    table = _filter_crosswalk_table(
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

    default_clids = merged_obs[data_label_column].map(create_temp_asctb_id)
    _set_defaults(merged_obs, data_clid_column, default_clids)
    _set_defaults(merged_obs, table_clid_label_column, merged_obs[data_label_column])
    _set_defaults(merged_obs, data_match_column, "skos:exactMatch")

    result = matrix.copy()
    result.obs = merged_obs
    _fix_obs_columns_dtype(result)

    return result


def _filter_invalid_rows(matrix: anndata.AnnData, column: str) -> anndata.AnnData:
    """Filter out rows where the obs column's value is NaN or the empty string ''.

    Args:
        matrix (anndata.AnnData): Matrix to filter
        column (str): Column in obs

    Returns:
        anndata.AnnData: Filtered subset matrix
    """
    obs_subset = matrix.obs[column]
    mask = obs_subset.notna() & (obs_subset != "")
    return matrix[mask, :]


def _filter_crosswalk_table(
    table: pd.DataFrame,
    organ_id: str,
    organ_level: str,
    organ_id_column: str,
    organ_level_column: str,
    table_label_column: str,
) -> pd.DataFrame:
    """Filter the crosswalk table to only include rows with organ id and level.

    Also removes empty rows.

    Args:
        table (pd.DataFrame): Original full crosswalk table

    Returns:
        pd.DataFrame: Filtered table
    """
    organ_id_rows = table[organ_id_column].str.lower() == organ_id.lower()
    organ_level_rows = table[organ_level_column].str.lower() == organ_level.lower()
    filtered_table = table[organ_id_rows & organ_level_rows]
    normalized_table = filtered_table.dropna(how="all")
    unique_table = normalized_table.drop_duplicates(table_label_column)
    return unique_table


def _set_defaults(
    obs: pd.DataFrame, column: str, defaults: t.Union[pd.Series, str]
) -> None:
    """Replace nan values with defaults in a column.

    Args:
        obs (pd.DataFrame): Data frame
        column (str): Column to update
        defaults (t.Union[pd.Series, str]): Default values
    """
    obs.loc[obs[column].isna(), column] = defaults


def _fix_obs_columns_dtype(matrix: anndata.AnnData):
    """Converts object and category columns to string to prevent errors when writing h5ad file.

    Args:
      matrix (AnnData): Matrix to update
    """
    for column in matrix.obs.columns:
        array = matrix.obs[column]
        if array.dtype in ("category", "object"):
            matrix.obs[column] = array.astype(str)


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
            args.crosswalk_table_clid_label_column,
            args.crosswalk_table_match_column,
        ]
    )


def _is_header_row(row: t.List[str], args: argparse.ArgumentParser) -> bool:
    """Tests whether a row is the header row.

    Args:
        row (t.List[str]): Row to test
        args (argparse.ArgumentParser): Same arguments as provided to `main`

    Returns:
        bool: True if it is the header row
    """
    to_lower_set = lambda items: set(map(str.lower, items))
    columns = [
        args.crosswalk_table_organ_id_column,
        args.crosswalk_table_organ_level_column,
        args.crosswalk_table_label_column,
        args.crosswalk_table_clid_column,
        args.crosswalk_table_clid_label_column,
        args.crosswalk_table_match_column,
    ]

    return to_lower_set(columns).issubset(to_lower_set(row))


def _read_table(args: argparse.Namespace) -> t.Optional[pd.DataFrame]:
    """Read a crosswalking table. Metadata rows before the header are skipped.

    Args:
        args (argparse.Namespace): Same arguments as provided to `main`

    Returns:
        pd.DataFrame: A data frame with the table data
    """
    with open(args.crosswalk_table) as file:
        for row in csv.reader(file):
            if _is_header_row(row, args):
                return pd.read_csv(file, names=row)
    return _get_empty_table(args)


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
        str(metadata["organ_id"]),
        str(metadata["organ_level"]),
        args.annotation_column,
        args.clid_column,
        args.match_column,
        _read_table(args),
        args.crosswalk_table_organ_id_column,
        args.crosswalk_table_organ_level_column,
        args.crosswalk_table_label_column,
        args.crosswalk_table_clid_column,
        args.crosswalk_table_clid_label_column,
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
        required=True,
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
        "--crosswalk-table-clid-label-column",
        default="CL_Label",
        help="Column with CL label in crosswalking table",
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
