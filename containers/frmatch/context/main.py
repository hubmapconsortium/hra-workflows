
import os
import sys
CODE_PATH = "../frmatch"
sys.path.insert(0, os.path.abspath(CODE_PATH))
import frmatch
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import pickle

import csv
import logging
import typing as t
from pathlib import Path

from src.algorithm import Algorithm, RunResult, add_common_arguments
from src.util.layers import set_data_layer

class FRmatchOrganMetadata(t.TypedDict):
    model: str

class FRmatchOptions(t.TypedDict):
    reference_data_dir: Path
    cluster_header: str
    gene_symbol_query: str
    ensembl_id_query: str

class FRmatchAlgorithm(Algorithm[FRmatchOrganMetadata, FRmatchOptions]):
    def __init__(self):
        super().__init__("predicted_labels")

    def do_run(
        self,
        file_h5ad: Path,
        organ: str,
        metadata: FRmatchOrganMetadata,
        options: FRmatchOptions,
    ) -> RunResult:
        """Annotate data using FRmatch."""
        adata_query = sc.read_h5ad(file_h5ad)

        cluster_header = options["cluster_header"]
        gene_symbol_query = options["gene_symbol_query"]
        ensembl_id_query = options["ensembl_id_query"]

        # Run leiden clustering
        sc.pp.neighbors(adata_query)
        sc.tl.leiden(adata_query, resolution = 0.1, key_added = cluster_header)
        print(adata_query.obs[cluster_header].value_counts())
        
        # Loading organ reference
        print("AAA")
        print(organ, options)
        adata_organ = self.get_reference(organ, options)

        print("BBB")
        # Subsetting adata_query to adata_organ feature space
        adata_query = self.unique_symbols(adata_query, gene_symbol_query, ensembl_id_query)
        print("CCC")
        adata_query = self.into_ref_space(adata_query, adata_organ)
        adata_query.var.index = adata_query.var[gene_symbol_query]

        print("running frmatch")
        # Run FRmatch cell2cluster
        settings, p_values = frmatch.FRmatch_cell2cluster(adata_query, adata_organ, subsamp_iter = 10, cluster_header_query = cluster_header, cluster_header_ref = cluster_header, save = True)
        # p_values = pd.read_pickle("frmatch_results_cell2cluster_querytoref.pkl")["results"]
        # run pval adj, force match but with lower confidence (unassigned)
        print("DDD")
        # Extracting annotation (per cell) from p_values
        print(p_values)
        annotation = self.get_annotation(p_values)
        print("EEE")
        print(annotation)
        for i in range(annotation.shape[0], adata_query.shape[0]): 
            temp = pd.DataFrame(dict(zip(annotation.columns, ["0", "index", "unassigned", 1.00])), index = [0])
            annotation = pd.concat([annotation, temp])
        annotation = annotation.rename(columns = {0: "frmatch_annotation", 1: "frmatch_confidence"})
        print(annotation)
        self.copy_annotations(adata_query, annotation)
        print(adata_query)
        return {"data": adata_query, 
        "organ_level": metadata["model"].replace(".", "_"), 
        "prediction_column": "frmatch_annotation"}

    def copy_annotations(
        self,
        matrix: ad.AnnData,
        annotation: pd.DataFrame,
    ) -> None:
        """Copies annotations from one matrix to another.

        Args:
            matrix (anndata.AnnData): Matrix to copy to
            annotated_matrix (anndata.AnnData): Matrix to copy from
        """
        matrix.obs = matrix.obs.merge(
            annotation,
            how="left",
            right_on="index",
            left_index=True,
            suffixes=(None, "_frmatch"),
        )
    def unique_symbols(self, query: ad.AnnData, gene_symbol_query:str, ensembl_id_query:str) -> ad.AnnData: 
        if ensembl_id_query == "index": 
            query.var["index"] = query.var.index
        query.var[ensembl_id_query + "_clean"] = [str(val).split(".")[0] for val in query.var[ensembl_id_query]]
        subset_query = query.var.copy().dropna(how = "any")
        if subset_query[subset_query.duplicated(ensembl_id_query + "_clean", keep = False)].shape[0] != 0: 
            print("WARNING: ensembl_id is not unique, need to debug")
        subset_query = subset_query.drop_duplicates(gene_symbol_query, keep = "first")
        query = query[:,list(subset_query.index)]
        query.var.index = query.var[gene_symbol_query]
        return query

    def get_reference(self, organ: str, options: FRmatchOptions) -> ad.AnnData: 
        print("entered get_reference")
        if organ == "lung" or organ == "UBERON:0002048": 
            organ = "hlca"
            print(options["reference_data_dir"])
            reference_data_path = self.find_reference_data(options["reference_data_dir"], organ)
            adata = sc.read_h5ad(reference_data_path)
            print(adata)
            return adata
    
    def find_reference_data(self, dir: Path, organ: str) -> Path:
        """Finds the reference data directory for an organ.

        Args:
            dir (Path): Directory to search
            organ (str): Organ name

        Raises:
            ValueError: If no reference data could be found

        Returns:
            Path: The data directory
        """
        print("entered find_reference_data")
        def is_reference_data_candidate(path: Path):
            print("entered is_reference_data_candidate")
            print(organ.lower)
            print(path.stem.lower())
            return (
                path.is_file()
                and path.suffix == ".h5ad"
                and organ.lower() in path.stem.lower()
            )

        return self._find_in_dir(
            dir,
            is_reference_data_candidate,
            f"Cannot find reference data for organ '{organ}'",
            f"Multiple reference data candidates for organ '{organ}'",
        )

    def _find_in_dir(
        self, dir: Path, cond: t.Callable[[Path], bool], error_msg: str, warn_msg: str
    ) -> Path:
        """Search a directory for a entry which passes the provided test.

        Args:
            dir (Path): Directory to search
            cond (t.Callable[[Path], bool]): Test used to match sub entries
            error_msg (str): Error message used when no entries match
            warn_msg (str): Warning message use when multiple entries match

        Raises:
            ValueError: If there are no matching sub entries

        Returns:
            Path:
                The matching entry.
                If multiple entries match the one with the shortest name is returned.
        """
        print("endered _find_in_dir")
        print(dir, cond, error_msg, warn_msg)
        print(os.getcwd())
        print(dir.iterdir())
        candidates = list(filter(cond, dir.iterdir()))
        print(candidates)
        candidates.sort(key=lambda path: len(path.name))
        print(candidates)
        if not candidates:
            raise ValueError(error_msg)
        elif len(candidates) > 1:
            warn(warn_msg)
        return candidates[0]

    def into_ref_space(self, query: ad.AnnData, ref: ad.AnnData) -> ad.AnnData: 
        genes_in_common = list(set(ref.var.index).intersection(set(query.var.index)))
        query = query[:, genes_in_common]
        return query

    def get_annotation(self, results: pd.DataFrame) -> pd.DataFrame: 
        print("entered get_annotation")
        p_adj_method = "fdr_by"
        sig_level = 0.05

        # Pivotting results to ref_cluster as index and query_cluster as columns
        results = results.sort_values(["query_cluster", "index", "p_value"], ascending = [True, True, False])
        results = results[~results.duplicated(["index", "query_cluster", "ref_cluster"], keep = "first")]
        results = results.pivot(index = "ref_cluster", columns = ["query_cluster", "index"], values = "p_value").replace(np.nan, 0)
        print("results", results.shape)
        # Adjusting p-values
        results_2 = frmatch.padj_FRmatch(results, p_adj_method = p_adj_method)
        print("results_2", results_2.shape)
        # Setting values to zero unless passes sig_level threshold
        results_3 = results_2.applymap(lambda x: x if x > sig_level else 0)
        print("results_3", results_3.shape)
        # I think we need to keep the unassigned row b/c that occurs when the sum of adjusted p-values = 0
        # Adding the unassigned row
        samples = results_3.columns[results_3.sum() == 0]
        temp = pd.DataFrame(dict(zip(samples, [1] * len(samples))), index = ["unassigned"])
        if temp.shape[1] != 0: 
            temp.columns.names = ["query_cluster", "index"]
            temp.index.name = "ref_cluster"
            results_4 = pd.concat([results_3, temp])
        else: 
            results_4 = results_3.copy()
        print("results_4", results_4.shape)
        # idxmax: returns ref_cluster with largest value per sample
        temp = pd.DataFrame(results_4.max())
        annotation = pd.DataFrame(results_4.idxmax()).reset_index()
        annotation[1] = list(temp[0])
#         annotation = zip(annotation["index"], annotation[0], annotation[1])
        # do we want to return a dataframe or zip?
        return annotation

def _get_arg_parser():
    parser = add_common_arguments()

    parser.add_argument(
        "--reference-data-dir",
        type=Path,
        required=True,
        help="Path to directory with reference data",
    )
    # parser.add_argument("--query-layers-key", help="Data layer to use")

    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    print(parser)
    args = parser.parse_args()
    algorithm = FRmatchAlgorithm()
    options = {'cluster_header': 'ann_finest_level', 'gene_symbol_query': 'hugo_symbol', 'ensembl_id_query': 'index'}
    result = algorithm.run(**args.__dict__, **options)
    result.save()
