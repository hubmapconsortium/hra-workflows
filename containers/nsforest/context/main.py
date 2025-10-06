import typing as t
import pandas as pd
import anndata
from nsforest import ns, nsforesting

from src.util.helper_gene_expression import get_arg_parser, common_main


def format_nsforest_markers_df(df: pd.DataFrame, cluster_header: str) -> pd.DataFrame:
    """Format the output of NSForest into a dataframe with cell type and marker genes as columns.
    This function mimics the exact structure of format_marker_genes_df from gene-expression.

    Args:
        df (pd.DataFrame): Output from NSForest
        cluster_header (str): Cell type id column

    Returns:
        pd.DataFrame: A data frame with `cluster_header` and marker_genes columns
    """

    # Ensure we have the right column names
    if "NSForest_markers" in df.columns:
        df = df.rename(columns={"NSForest_markers": "marker_genes"})

    # Convert marker_genes to lists if they aren't already
    df["marker_genes"] = df["marker_genes"].apply(
        lambda x: x if isinstance(x, list) else [x] if pd.notna(x) else []
    )

    # Keep the original cluster_header values
    df = df[[cluster_header, "marker_genes"]].copy()

    return df


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

    try:
        # NSForest preprocessing
        matrix = ns.pp.prep_medians(matrix, cluster_header)
        matrix = ns.pp.prep_binary_scores(matrix, cluster_header)

        markers = nsforesting.NSForest(
            matrix,
            cluster_header,
            save=False,
            save_supplementary=False,
            output_folder=".",
            outputfilename_prefix="nsforest",
        )

        if isinstance(markers, pd.DataFrame):
            df = markers.copy()
            df = df.rename(columns={"clusterName": cluster_header})
        else:
            df = pd.DataFrame(markers)
            df = df.rename(columns={"clusterName": cluster_header})

        # Format into markers dataframe using the same structure as gene-expression
        markers_df = format_nsforest_markers_df(df, cluster_header)

        return markers_df[[cluster_header, "marker_genes"]]
    
    except Exception as e:
        print(f"ERROR: NSForest failed with error: {e}")
        empty_df = pd.DataFrame(columns=[cluster_header, "marker_genes"])
        empty_df["marker_genes"] = empty_df["marker_genes"].astype(object)
        return empty_df


def main():
    """Main function for NSForest marker gene analysis."""
    parser = get_arg_parser("Add NSForest markers to h5ad data", "get_marker_nsforest")
    args = parser.parse_args()

    # Use the common main function with NSForest-specific marker function
    common_main(args, get_marker_nsforest)


if __name__ == "__main__":
    main()
