import typing as t
import pandas as pd
import anndata
from nsforest import ns, nsforesting
from pathlib import Path

from src.util.helper_gene_expression import get_arg_parser, common_main


def get_marker_nsforest(
    matrix: anndata.AnnData, cluster_header: str, n_genes: int
) -> pd.DataFrame:
    """Run NSForest and return a DataFrame with **cluster_header** and **marker_genes**.

    Args:
        matrix (anndata.AnnData): Filtered matrix with ensemble data
        cluster_header (str): Column name for cell type grouping
        n_genes (int): Number of top genes per cell type (ignored by NSForest)

    Returns:
        pd.DataFrame: DataFrame with columns [cluster_header, marker_genes]
    """
    adata = matrix.copy()

    # NSForest preprocessing
    adata = ns.pp.prep_medians(adata, cluster_header)
    adata = ns.pp.prep_binary_scores(adata, cluster_header)

    markers = nsforesting.NSForest(
        adata,
        cluster_header,
        save=False,
        save_supplementary=False,
        output_folder=".",
        outputfilename_prefix="nsforest",
    )

    if isinstance(markers, pd.DataFrame):
        df = markers.copy()
        df.rename(columns={"clusterName": cluster_header}, inplace=True)
    else:
        df = pd.DataFrame(markers)
        df.rename(columns={"clusterName": cluster_header}, inplace=True)

    if "NSForest_markers" in df.columns:
        df.rename(columns={"NSForest_markers": "marker_genes"}, inplace=True)

    return df[[cluster_header, "marker_genes"]]


def main():
    """Main function for NSForest marker gene analysis."""
    parser = get_arg_parser("Add NSForest markers to h5ad data")
    args = parser.parse_args()

    # Use the common main function with NSForest-specific marker function
    common_main(args, get_marker_nsforest)


if __name__ == "__main__":
    main()
