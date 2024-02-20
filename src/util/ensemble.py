import os
import re
import typing as t

import anndata
import pandas as pd

from .ids import create_gene_id

_ENSEMBLE_COLUMN = "ensemble"
_GENE_COLUMN = "hgnc"
_GENE_NAME_COLUMN = "gene_name"

_ENSEMBLE_VERSION_REPLACE_REGEX = re.compile(r"\.\d+$")


def strip_version(value: str) -> str:
    """Removes the version specifier from an ensembl id.

    Args:
        value (str): Id with optional version specified

    Returns:
        str: Ensembl id without version
    """
    return re.sub(_ENSEMBLE_VERSION_REPLACE_REGEX, "", value)


def add_ensemble_data(
    matrix: anndata.AnnData, ensemble: t.Union[str, bytes, os.PathLike, pd.DataFrame]
) -> anndata.AnnData:
    """Add ensembl and gene information to `var` from a lookup file.

    Args:
        matrix (anndata.AnnData): Original matrix
        ensemble (t.Union[str, bytes, os.PathLike, pd.DataFrame]): Ensemble lookup or path to csv file

    Returns:
        anndata.AnnData: Matrix with ensembl information
    """
    if not isinstance(ensemble, pd.DataFrame):
        ensemble = pd.read_csv(ensemble, dtype=str)
    ensemble = ensemble.drop_duplicates(_ENSEMBLE_COLUMN)

    index = matrix.var.index
    keys = index.map(strip_version)
    merged_var = matrix.var.merge(
        ensemble, how="left", left_on=keys, right_on=_ENSEMBLE_COLUMN
    )
    merged_var.index = index
    merged_var[_GENE_COLUMN].fillna(index.to_series().map(create_gene_id), inplace=True)
    merged_var[_GENE_NAME_COLUMN].fillna(index.to_series(), inplace=True)

    result = matrix.copy()
    result.var = merged_var
    return result
