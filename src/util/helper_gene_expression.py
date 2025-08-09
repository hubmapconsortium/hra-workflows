import argparse
import json
import traceback
import typing as t
from pathlib import Path
import anndata
import numpy as np
import pandas as pd
import scanpy as sc

from src.util.ensemble import add_ensemble_data


MIN_CELLS_PER_CT = 2  # need a filter to remove CTs with one cell as sc.tl.rank_genes_groups gives an error


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
    
    # Get cells of this type
    cells_of_type = matrix.obs[matrix.obs[clid_column] == cell_type]

    cell_indices = [
        matrix.obs.index.get_loc(cell_index)
        for cell_index in cells_of_type.index
    ]
    gene_loc = matrix.var.index.get_loc(gene)
    expr_values = matrix.X[cell_indices, gene_loc]
    mean_expr = expr_values.mean()

    return float(mean_expr)


def get_marker_genes_with_expr(
    matrix: anndata.AnnData,
    gene_data: pd.DataFrame,
    clid_column: str,
    cell_type: str,
    marker_genes: t.Iterable[str],
) -> t.List[dict]:
    """Get the mean expression for all marker genes.

    Args:
        matrix (anndata.AnnData): Data
        gene_data (pd.DataFrame): Data frame with additional gene data indexed by gene
        clid_column (str): Cell type id column
        cell_type (str): Cell type name
        marker_genes: Iterable of gene names

    Returns:
        t.List[dict]: Mean expressions, the dicts contains "gene_label" and "mean_gene_expr_value"
    """
    
    output: list[dict] = []

    for i, gene in enumerate(marker_genes):
        try:
            # Check if gene exists in gene_data
            if gene in gene_data.index:
                data = gene_data.loc[gene]
                output.append(
                    {
                        "gene_id": data.get("hgnc", ""),
                        "gene_label": data.get("gene_name", gene),
                        "ensembl_id": gene,
                        "mean_gene_expr_value": get_mean_expr_value(
                            matrix, clid_column, cell_type, gene
                        ),
                    }
                )
            else:
                # Gene not found in gene_data, use default values
                output.append(
                    {
                        "gene_id": "",
                        "gene_label": gene,
                        "ensembl_id": gene,
                        "mean_gene_expr_value": get_mean_expr_value(
                            matrix, clid_column, cell_type, gene
                        ),
                    }
                )
        except Exception as e:
            print(f"ERROR: Exception type: {type(e)}")
            continue
    
    return output


def prepare_filtered_matrix(
    matrix: anndata.AnnData, clid_column: str, ensemble: Path, min_cell_types: int = 2
) -> anndata.AnnData:
    """Prepare filtered matrix with ensemble data for marker gene analysis.

    Args:
        matrix (anndata.AnnData): Original matrix
        clid_column (str): Cell type id column
        ensemble (Path): Path to ensemble lookup data
        min_cell_types (int): Minimum number of unique cell types required

    Returns:
        anndata.AnnData: Filtered matrix with ensemble data

    Raises:
        ValueError: If too few unique cell types for marker discovery
    """
    matrix.raw = matrix  

    filtered_matrix = filter_matrix(matrix, clid_column)
    filtered_matrix = add_ensemble_data(filtered_matrix, ensemble)

    if len(filtered_matrix.obs[clid_column].unique()) <= min_cell_types - 1:
        raise ValueError(
            f"Too few unique cell types for marker discovery (found {len(filtered_matrix.obs[clid_column].unique())}, need at least {min_cell_types})"
        )

    return filtered_matrix


def merge_markers_into_matrix(
    matrix: anndata.AnnData,
    markers_df: pd.DataFrame,
    clid_column: str,
    gene_expr_column: str,
) -> anndata.AnnData:
    """Merge marker gene data into the matrix observations.

    Args:
        matrix (anndata.AnnData): Original matrix
        markers_df (pd.DataFrame): DataFrame with marker genes and expression data
        clid_column (str): Cell type id column name in matrix.obs
        gene_expr_column (str): Column name for gene expression data

    Returns:
        anndata.AnnData: Updated matrix with marker gene data
    """
    
    # Merge into obs - explicitly specify the merge column to avoid unhashable type errors
    merged_obs = matrix.obs.merge(
        markers_df[[clid_column, gene_expr_column]],
        how="left",
        left_on=clid_column,
        right_on=clid_column,
    )

    
    merged_obs.fillna({gene_expr_column: "[]"}, inplace=True)

    
    # Handle the case where pandas merge created suffixes
    actual_gene_expr_column = gene_expr_column
    if gene_expr_column not in merged_obs.columns:
        # Check for suffixed columns
        if f"{gene_expr_column}_x" in merged_obs.columns and f"{gene_expr_column}_y" in merged_obs.columns:
            # Use the _y column (from markers_df) and drop the _x column (from matrix.obs)
            actual_gene_expr_column = f"{gene_expr_column}_y"
            merged_obs = merged_obs.drop(columns=[f"{gene_expr_column}_x"])
        elif f"{gene_expr_column}_y" in merged_obs.columns:
            actual_gene_expr_column = f"{gene_expr_column}_y"
        else:
            raise KeyError(f"Column {gene_expr_column} not found in merged_obs")
    
    
    if actual_gene_expr_column in merged_obs.columns:
        
        try:
            merged_obs[actual_gene_expr_column] = merged_obs[actual_gene_expr_column].apply(json.dumps)

        except Exception as e:
            # Handle the case where json.dumps fails
            merged_obs[actual_gene_expr_column] = merged_obs[actual_gene_expr_column].apply(lambda x: str(x) if x is not None else "[]")




        # Rename back to original column name if needed
        if actual_gene_expr_column != gene_expr_column:
            merged_obs = merged_obs.rename(columns={actual_gene_expr_column: gene_expr_column})

    else:
        raise KeyError(f"Column {actual_gene_expr_column} not found in merged_obs")
    
    merged_obs.index = matrix.obs.index
    
    matrix.obs = merged_obs
    
    return matrix


def get_gene_expr_with_marker_function(
    matrix: anndata.AnnData,
    ensemble: Path,
    clid_column: str,
    gene_expr_column: str,
    n_genes: int,
    marker_function: t.Callable[[anndata.AnnData, str, int], pd.DataFrame],
) -> anndata.AnnData:
    """Generic function to compute and add gene mean expressions for all cells.

    Args:
        matrix (anndata.AnnData): Data
        ensemble (Path): Ensemble lookup data path
        clid_column (str): Cell type id column
        gene_expr_column (str): Column to store gene expression on
        n_genes (int): Number of top genes per cell type to save
        marker_function (callable): Function that takes (filtered_matrix, clid_column, n_genes) and returns markers_df

    Returns:
        anndata.AnnData: Updated data with gene expressions
    """

    filtered_matrix = prepare_filtered_matrix(matrix, clid_column, ensemble)
    
    markers_df = marker_function(filtered_matrix, clid_column, n_genes)

    # Compute mean expression per marker gene
    
    # Handle empty DataFrame case
    if len(markers_df) == 0:
        markers_df[gene_expr_column] = []
    else:
        markers_df[gene_expr_column] = markers_df.apply(
            lambda row: get_marker_genes_with_expr(
                matrix,
                filtered_matrix.var,
                clid_column,
                row[clid_column],
                row["marker_genes"],
            ),
            axis=1,
        )
    print(f"Mean expressions computed successfully")

    result = merge_markers_into_matrix(
        matrix,
        markers_df,
        clid_column,
        gene_expr_column,
    )

    
    return result


def get_arg_parser(description: str, marker_function_name: str = None) -> argparse.ArgumentParser:
    """Common argument parser for both implementations."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("matrix", type=anndata.read_h5ad, help="h5ad data file")
    parser.add_argument(
        "--ensemble-lookup",
        default="/src/assets/ensemble-lookup.csv",
        help="Ensemble lookup data csv",
    )
    parser.add_argument("--clid-column", default="clid", help="Column with cell id")
    parser.add_argument(
        "--gene-expr-column", default="gene_expr", help="Column for gene_expr"
    )
    parser.add_argument(
        "--gene-expr-count",
        type=int,
        default=200,
        help="number of top genes per cell type to save",
    )
    
    # Determine output filename based on marker function name
    if marker_function_name:
        if "nsforest" in marker_function_name.lower():
            default_output = "matrix_with_nsforest.h5ad"
        elif "scanpy" in marker_function_name.lower() or "gene" in marker_function_name.lower():
            default_output = "matrix_with_gene_expr.h5ad"
        else:
            default_output = "matrix_with_gene_expr.h5ad"
    else:
        default_output = "matrix_with_gene_expr.h5ad"
    
    parser.add_argument(
        "--output-matrix",
        type=Path,
        default=default_output,
        help="matrix with gene expressions output path",
    )
    parser.add_argument(
        "--output-report",
        type=Path,
        default="report.json",
        help="Report file path in case of errors",
    )
    return parser


def common_main(
    args: argparse.Namespace,
    marker_function: t.Callable[[anndata.AnnData, str, int], pd.DataFrame],
) -> None:
    """Common main function for both implementations."""
    try:
        matrix = get_gene_expr_with_marker_function(
            args.matrix,
            args.ensemble_lookup,
            args.clid_column,
            args.gene_expr_column,
            args.gene_expr_count,
            marker_function, 
        )

        matrix.write_h5ad(args.output_matrix)
        
        # Create success report
        args.output_report.parent.mkdir(parents=True, exist_ok=True)
        with args.output_report.open("w") as fh:
            json.dump(
                {
                    "status": "success",
                },
                fh,
                indent=4,
            )

    except Exception as error:
        print("Exception occurred:", error)
        args.output_report.parent.mkdir(parents=True, exist_ok=True)
        with args.output_report.open("w") as fh:
            json.dump(
                {
                    "status": "failure",
                    "cause": repr(error),
                    "traceback": traceback.format_tb(error.__traceback__),
                },
                fh,
                indent=4,
            )
