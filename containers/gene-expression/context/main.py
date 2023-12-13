import argparse
import json
import typing as t
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import scanpy as sc

MIN_CELLS_PER_CT = 2  # need a filter to reomve CTs with one cell as sc.tl.rank_genes_groups gives an error


def filter_matrix(matrix: anndata.AnnData, clid_column: str) -> anndata.AnnData:
    """Filters data to only include cells and cell types with at least `MIN_CELLS_PER_CT`.

    Args:
        matrix (anndata.AnnData): Data to filter
        clid_column (str): Cell type id column

    Returns:
        anndata.AnnData: Filtered data
    """
    ct_counts = matrix.obs[clid_column].value_counts()
    valid_cts = ct_counts.index[ct_counts >= MIN_CELLS_PER_CT]
    mask = np.isin(matrix.obs[clid_column], valid_cts)
    return matrix[mask, :]


def format_marker_genes_df(df: pd.DataFrame, clid_column: str) -> pd.DataFrame:
    """Format the output of `scanpy.tl.rank_genes_groups` into a dataframe with
    cell type and marker genes as columns.

    Args:
        df (pd.DataFrame): Output from `scanpy.tl.rank_genes_groups`
        clid_column (str): Cell type id column

    Returns:
        pd.DataFrame: A data frame with `clid_column` and marker_genes columns
    """
    df = df.transpose()
    df["marker_genes"] = df.apply(lambda row: row.tolist(), axis=1)
    df = df["marker_genes"].rename_axis(clid_column).reset_index()
    return df


def get_mean_expr_value(
    matrix: anndata.AnnData, clid_column: str, cell_type: str, gene: str
) -> float:
    """Computes the mean expression for a cell type and gene.

    Args:
        matrix (anndata.AnnData): Data
        clid_column (str): Cell type id column
        cell_type (str): Cell type name
        gene (str): Gene name

    Returns:
        float: The mean expression
    """
    cell_indices = [
        matrix.obs.index.get_loc(cell_index)
        for cell_index in matrix.obs[matrix.obs[clid_column] == cell_type].index
    ]
    mean_expr = matrix.X[cell_indices, matrix.var.index.get_loc(gene)].mean()
    return float(mean_expr)


def get_marker_genes_with_expr(
    matrix: anndata.AnnData, clid_column: str, cell_type: str, marker_genes: str
) -> t.List[dict]:
    """Get the mean expression for all marker genes.

    Args:
        matrix (anndata.AnnData): Data
        clid_column (str): Cell type id column
        cell_type (str): Cell type name
        gene (str): Gene name

    Returns:
        t.List[dict]: Mean expressions, the dicts contains "gene_label" and "mean_gene_expr_value"
    """
    output = []
    for gene in marker_genes:
        output.append(
            {
                "gene_label": gene,
                "mean_gene_expr_value": get_mean_expr_value(
                    matrix, clid_column, cell_type, gene
                ),
            }
        )
    return output


def get_gene_expr(
    matrix: anndata.AnnData, clid_column: str, gene_expr_column: str
) -> anndata.AnnData:
    """Computes and adds gene mean expressions for all cells in the annotated data.

    Args:
        matrix (anndata.AnnData): Data
        clid_column (str): Cell type id column
        gene_expr_column (str): Column to store gene expression on

    Returns:
        anndata.AnnData: Updated data with gene expressions
    """
    matrix.raw = matrix  # for getting gene names as output for sc.tl.rank_genes_groups
    filtered_matrix = filter_matrix(matrix, clid_column)
    sc.tl.rank_genes_groups(filtered_matrix, groupby=clid_column, n_genes=10)
    ct_marker_genes_df = format_marker_genes_df(
        pd.DataFrame(filtered_matrix.uns["rank_genes_groups"]["names"]),
        clid_column,
    )

    ct_marker_genes_df[gene_expr_column] = ct_marker_genes_df.apply(
        lambda row: get_marker_genes_with_expr(
            matrix, clid_column, row[clid_column], row["marker_genes"]
        ),
        axis=1,
    )
    merged_obs = matrix.obs.merge(
        ct_marker_genes_df[[clid_column, gene_expr_column]], how="left"
    )
    merged_obs.fillna({gene_expr_column: "[]"}, inplace=True)
    merged_obs[gene_expr_column] = merged_obs[gene_expr_column].apply(json.dumps)
    merged_obs.index = matrix.obs.index
    matrix.obs = merged_obs
    return matrix


def main(args: argparse.Namespace):
    """Computes gene mean expression for all cells in an annotated h5ad file and
    saves the result to another h5ad file.

    Args:
        args (argparse.Namespace):
            CLI arguments, must contain "matrix", "clid_column",
            "gene_expr_column", and "output_matrix"
    """
    matrix = get_gene_expr(args.matrix, args.clid_column, args.gene_expr_column)
    matrix.write_h5ad(args.output_matrix)


def _get_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Add gene expressions to h5ad data")
    parser.add_argument("matrix", type=anndata.read_h5ad, help="h5ad data file")
    parser.add_argument("--clid-column", default="clid", help="Column with cell id")
    parser.add_argument(
        "--gene-expr-column", default="gene_expr", help="Column for gene_expr"
    )
    parser.add_argument(
        "--output-matrix",
        type=Path,
        default="matrix_with_gene_expr.h5ad",
        help="matrix with gene expressions output path",
    )

    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
