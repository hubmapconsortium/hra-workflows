import typing as t
import pandas as pd
import anndata
import scanpy as sc
from pathlib import Path

from src.util.helper_gene_expression import get_arg_parser, common_main


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


def get_marker_gene_scanpy(
    matrix: anndata.AnnData, cluster_header: str, n_genes: int
) -> pd.DataFrame:
    """Run scanpy rank_genes_groups and return a DataFrame with cluster and marker_genes columns.

    Args:
        matrix (anndata.AnnData): Filtered matrix with ensemble data
        cluster_header (str): Column name for cell type grouping
        n_genes (int): Number of top genes per cell type

    Returns:
        pd.DataFrame: DataFrame with columns [cluster_header, marker_genes]
    """
    adata = matrix

    # Run scanpy rank genes groups
    sc.tl.rank_genes_groups(adata, groupby=cluster_header, n_genes=n_genes)

    # De-fragment the rank_genes_groups results
    for key in ["names", "scores", "pvals", "pvals_adj", "logfoldchanges"]:
        if key in adata.uns["rank_genes_groups"]:
            adata.uns["rank_genes_groups"][key] = adata.uns["rank_genes_groups"][
                key
            ].copy()

    # Format into markers dataframe
    markers_df = format_marker_genes_df(
        pd.DataFrame(adata.uns["rank_genes_groups"]["names"]),
        cluster_header,
    )

    return markers_df[[cluster_header, "marker_genes"]]


def main():
    """Main function for gene expression marker gene analysis."""
    parser = get_arg_parser("Generate gene expression JSON from marker genes")
    args = parser.parse_args()

    # Use the common main function with gene expression-specific marker function
    common_main(args, get_marker_gene_scanpy)


if __name__ == "__main__":
    main()
