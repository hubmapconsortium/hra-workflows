import os
import re
import typing as t

import anndata
import pandas as pd


_ENSEMBLE_COLUMN = "ensemble"
_GENE_COLUMN = "hgnc"
_GENE_NAME_COLUMN = "gene_name"

_ENSEMBLE_VERSION_REPLACE_REGEX = re.compile(r"\.\d+$")


def strip_version(value: str) -> str:
    return re.sub(_ENSEMBLE_VERSION_REPLACE_REGEX, "", value)


def add_ensemble_data(
    matrix: anndata.AnnData, ensemble: t.Union[str, bytes, os.PathLike, pd.DataFrame]
) -> anndata.AnnData:
    if not isinstance(ensemble, pd.DataFrame):
        ensemble = pd.read_csv(ensemble, dtype=str)
    ensemble = ensemble.drop_duplicates(_ENSEMBLE_COLUMN)

    index = matrix.var.index
    keys = index.map(strip_version)
    merged_var = matrix.var.merge(
        ensemble, how="left", left_on=keys, right_on=_ENSEMBLE_COLUMN
    )
    merged_var.index = index
    merged_var[_GENE_COLUMN].fillna(index.to_series().map(_create_default_gene_id), inplace=True)
    merged_var[_GENE_NAME_COLUMN].fillna(index.to_series(), inplace=True)

    result = matrix.copy()
    result.var = merged_var
    return result


def _create_temp_asctb_id(value: str) -> str:
    """Create a temporary IRI based on a label.

    Args:
        value (str): Label for the row

    Returns:
        str: Temporary IRI
    """
    suffix = value.lower().strip()
    suffix = re.sub(r"\W+", "-", suffix)
    suffix = re.sub(r"[^a-z0-9-]+", "", suffix)
    return "ASCTB-TEMP:" + suffix


def _create_default_gene_id(value: str) -> str:
    if value.lower().startswith('ens'):
        return 'ensembl:' + value
    return _create_temp_asctb_id(value)
