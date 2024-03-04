import typing as t

import anndata


def set_data_layer(matrix: anndata.AnnData, layer: t.Optional[str]) -> anndata.AnnData:
    """Sets the active data layer.
    If the layer does not exist it is ignored.

    Args:
        matrix (anndata.AnnData): Original matrix with layers
        layer (t.Optional[str]): Name of layer or 'X' or 'raw'

    Returns:
        anndata.AnnData: A new matrix with the layer set as the X data matrix
    """
    if layer in ('X', None):
        return matrix
    
    matrix = matrix.copy()
    if layer == 'raw' and matrix.raw is not None:
        matrix.X = matrix.raw.X
    elif layer in matrix.layers:
        matrix.X = matrix.layers[layer].copy()

    return matrix
