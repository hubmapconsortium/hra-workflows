# TODO Dockerize!!!!!

# from calendar import c
import anndata
import numpy as np
import scanpy as sc
import pandas as pd
from algorithm import add_common_arguments
from pathlib import Path

#need a filter to reomve CTs with one cell as sc.tl.rank_genes_groups gives an error
MIN_CELLS_PER_CT = 2
CELL_TYPE_COLUMN = "predicted.ann_finest_level" #from annotations.csv
COUNT_COLUMN = "count" #from annotations.csv

def get_marker_genes(output_matrix: Path, output_annotations: Path):
    adata = anndata.read_h5ad(output_matrix)
    ct_counts = pd.read_csv(output_annotations)
    # adata.raw = adata

    valid_cts = ct_counts.loc[ct_counts[COUNT_COLUMN] >= MIN_CELLS_PER_CT,CELL_TYPE_COLUMN]
    mask = np.isin(adata.obs[CELL_TYPE_COLUMN], valid_cts)
    filtered_adata = adata[mask, :]

    sc.tl.rank_genes_groups(filtered_adata, groupby=CELL_TYPE_COLUMN, n_genes=10)

    ct_marker_genes = pd.DataFrame(filtered_adata.uns['rank_genes_groups']['names']).transpose()
    ct_marker_genes['marker_genes'] = ct_marker_genes.apply(lambda row: row.tolist(), axis=1)
    ct_marker_genes = ct_marker_genes['marker_genes'].rename_axis(CELL_TYPE_COLUMN).reset_index()

    def get_expr_value(cell_type,marker_genes):
        output = []
        cell_indices = [filtered_adata.obs.index.get_loc(cell_index) for cell_index in filtered_adata.obs[filtered_adata.obs[CELL_TYPE_COLUMN] == ct].index]
        for gene in marker_genes:
            output.append({'gene_label':gene,'mean_gene_expr_value':filtered_adata.X[cell_indices,filtered_adata.var.index.get_loc(gene)].mean()})
        return output

    ct_marker_genes['gene_expr'] = ct_marker_genes.apply(lambda row: get_expr_value(row['ct'], row['marker_genes']), axis=1)

    return ct_counts.merge(ct_marker_genes[[CELL_TYPE_COLUMN,'gene_expr']],how='left')


if __name__ == "__main__":
    parser = add_common_arguments()
    args = parser.parse_args()
    result = get_marker_genes(**args.__dict__)
    result.save()
