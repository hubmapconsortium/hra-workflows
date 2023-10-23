import argparse
import anndata
import numpy as np
import scanpy as sc
import pandas as pd
from pathlib import Path

#need a filter to reomve CTs with one cell as sc.tl.rank_genes_groups gives an error
MIN_CELLS_PER_CT = 2

def get_marker_genes(matrix: anndata.AnnData, annotation_column: str, gene_expr_column: str):
    #for getting gene names as output for sc.tl.rank_genes_groups
    matrix.raw = matrix

    ct_counts = matrix.obs[annotation_column].value_counts()
    valid_cts = ct_counts.index[ct_counts >= MIN_CELLS_PER_CT]
    mask = np.isin(matrix.obs[annotation_column], valid_cts)
    filtered_matrix = matrix[mask, :]

    sc.tl.rank_genes_groups(filtered_matrix, groupby=annotation_column, n_genes=10)

    ct_marker_genes = pd.DataFrame(filtered_matrix.uns['rank_genes_groups']['names']).transpose()
    ct_marker_genes['marker_genes'] = ct_marker_genes.apply(lambda row: row.tolist(), axis=1)
    ct_marker_genes = ct_marker_genes['marker_genes'].rename_axis(annotation_column).reset_index()

    def get_expr_value(cell_type,marker_genes):
        output = []
        cell_indices = [filtered_matrix.obs.index.get_loc(cell_index) for cell_index in filtered_matrix.obs[filtered_matrix.obs[annotation_column] == cell_type].index]
        for gene in marker_genes:
            output.append({'gene_label':gene,'mean_gene_expr_value':filtered_matrix.X[cell_indices,filtered_matrix.var.index.get_loc(gene)].mean()})
        return output

    ct_marker_genes[gene_expr_column] = ct_marker_genes.apply(lambda row: get_expr_value(row[annotation_column], row['marker_genes']), axis=1)

    merged_obs = matrix.obs.merge(ct_marker_genes[[annotation_column,gene_expr_column]], how='left')
    merged_obs[gene_expr_column] = merged_obs[gene_expr_column].astype(str)
    merged_obs.index = matrix.obs.index
    matrix.obs = merged_obs

    return matrix

def main(args: argparse.Namespace):
    matrix = get_marker_genes(args.matrix, args.annotation_column, args.gene_expr_column)
    matrix.write_h5ad(args.output_matrix)

def _get_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Add gene expressions to h5ad data")
    parser.add_argument("matrix", type=anndata.read_h5ad, help="h5ad data file")
    parser.add_argument("--annotation-column", required=True, help="Column with annotations")
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