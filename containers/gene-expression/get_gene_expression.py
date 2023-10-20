# TODO Dockerize!!!!!

from calendar import c
import anndata
import numpy as np
import scanpy as sc
import pandas as pd

adata = anndata.read_h5ad('result.h5ad')
adata.raw = adata

#need a filter to reomve CTs with one cell as sc.tl.rank_genes_groups gives an error
MIN_CELLS_PER_CT = 2
CELL_TYPE_COLUMN = "predicted.ann_finest_level"

ct_counts = adata.obs[CELL_TYPE_COLUMN].value_counts()
valid_cts = ct_counts.index[ct_counts >= MIN_CELLS_PER_CT]
mask = np.isin(adata.obs[CELL_TYPE_COLUMN], valid_cts)

filtered_adata = adata[mask, :]

sc.tl.rank_genes_groups(filtered_adata, groupby=CELL_TYPE_COLUMN, n_genes=10)

ct_gene_expr = pd.DataFrame(filtered_adata.uns['rank_genes_groups']['names']).transpose()

ct_gene_expr['gene_expr'] = ct_gene_expr.apply(lambda row: row.tolist(), axis=1)

ct_gene_expr = ct_gene_expr['gene_expr'].rename_axis('ct').reset_index()

def get_gene_expr(ct,genes):
    output = []
    cell_indices = [filtered_adata.obs.index.get_loc(cell_index) for cell_index in filtered_adata.obs[filtered_adata.obs[CELL_TYPE_COLUMN] == ct].index]
    for gene in genes:
        output.append({'gene_label':gene,'mean_gene_expr_value':filtered_adata.X[cell_indices,filtered_adata.var.index.get_loc(gene)].mean()})
    return output

ct_gene_expr['gene_expr'] = ct_gene_expr.apply(lambda row: get_gene_expr(row['ct'], row['gene_expr']), axis=1)

ct_gene_expr.to_csv('gene_expr.csv')
