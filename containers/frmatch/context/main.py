import os
import sys
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
from src.util.ensemble import add_ensemble_data


class FRmatchOrganMetadata(t.TypedDict):
    model: str
    organ_level: str
    cluster_header: str


class FRmatchOptions(t.TypedDict):
    reference_data_dir: Path
    query_layers_key: t.Optional[str]
    ensemble_lookup: Path


class FRmatchAlgorithm(Algorithm[FRmatchOrganMetadata, FRmatchOptions]):
    def __init__(self):
        super().__init__("frmatch_annotation")

    def do_run(
        self,
        file_h5ad: Path,
        organ: str,
        metadata: FRmatchOrganMetadata,
        options: FRmatchOptions,
    ) -> RunResult:
        """Annotate data using FRmatch."""
        adata_query = sc.read_h5ad(file_h5ad)
        adata_query = set_data_layer(adata_query, options["query_layers_key"])

        # Getting .var columns
        adata_query = self.get_var_cols(adata_query, options["ensemble_lookup"])
        cluster_header = metadata["cluster_header"]

        # Loading organ reference
        adata_organ = self.get_reference(options, metadata)

        # Run leiden clustering
        sc.pp.neighbors(adata_query)
        sc.tl.leiden(adata_query, resolution=0.1, key_added=cluster_header)

        # Subsetting adata_query to adata_organ feature space
        adata_query = self.unique_symbols(adata_query.copy())
        adata_query = self.into_ref_space(adata_query.copy(), adata_organ.copy())

        # Run FRmatch cell2cluster
        settings, p_values = frmatch.FRmatch_cell2cluster(
            adata_query,
            adata_organ,
            subsamp_iter=2000,
            cluster_header_query=cluster_header,
            cluster_header_ref=cluster_header,
            save=True,
        )

        # Extracting annotation (per cell) from p_values
        annotation = self.get_annotation(p_values)
        annotation = annotation.rename(
            columns={0: "frmatch_annotation", 1: "frmatch_confidence"}
        )

        # Adding annotation to adata_query.obs
        adata_query = self.copy_annotations(adata_query, annotation)

        # Remove unlabeled cells
        adata_query.obs.index = adata_query.obs.index.astype(str)
        adata_query = adata_query[~adata_query.obs[self.prediction_column].astype(str).isin([False, "False", "unassigned"])]

        return {
            "data": adata_query,
            "model": metadata["model"],
            "organ_level": metadata["organ_level"],
            "cluster_header": metadata["cluster_header"],
            "prediction_column": "frmatch_annotation",
            "prediction_confidence": "frmatch_confidence",
        }

    def get_var_cols(self, adata: ad.AnnData, ensemble_lookup: Path) -> ad.AnnData:
        """Finding feature_name and ensembl_id columns in adata.var.

        Args:
            query (ad.AnnData): AnnData to rename .var columns.
            ensemble_lookup (Path): Path to ensemble lookup CSV file.

        Returns:
            ad.AnnData: AnnData with "feature_name" and "ensembl_id" in .var.
        """
        # Checking for feature_name
        feature_name = False
        for col in adata.var.columns:
            if "feature_name" in col or "hugo_symbol" in col or "gene_symbol" in col:
                feature_name = True
                break
        if feature_name:
            adata.var = adata.var.rename(columns={col: "feature_name"})
        else:
            logging.warning(f"Missing feature_name in adata.var column.")

        # Checking for ensembl_id
        ensembl_id = False
        for col in adata.var.columns:
            if "ensembl_id" in col:
                ensembl_id = True
                break
        if ensembl_id:
            adata.var = adata.var.rename(columns={col: "ensembl_id"})
        else:  # if 'ensembl_id' not in adata.var, must be unlabeled in index
            adata.var["ensembl_id"] = adata.var.index.astype(str)
        print("Before renaming")
        print(adata.var.columns)
        if not feature_name:
            adata = add_ensemble_data(adata, ensemble_lookup)
            adata.var = adata.var.rename(columns={"gene_name": "feature_name"})
        print("After renaming")
        print(adata.var.columns)
        return adata

    def get_reference(
        self, options: FRmatchOptions, metadata: FRmatchOrganMetadata
    ) -> ad.AnnData:
        reference_data_path = self.find_reference_data(
            options["reference_data_dir"], metadata["model"]
        )
        adata = sc.read_h5ad(reference_data_path)
        return adata

    def find_reference_data(self, dir: Path, model: str) -> Path:
        """Finds the reference data directory for a model.

        Args:
            dir (Path): Directory to search
            model (str): Organ name

        Raises:
            ValueError: If no reference data could be found

        Returns:
            Path: The data directory
        """

        def is_reference_data_candidate(path: Path):
            return (
                path.is_file()
                and path.suffix == ".h5ad"
                and model.lower() in path.stem.lower()
            )

        return self._find_in_dir(
            dir,
            is_reference_data_candidate,
            f"Cannot find reference data for organ '{model}'",
            f"Multiple reference data candidates for organ '{model}'",
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
        candidates = list(filter(cond, dir.iterdir()))
        candidates.sort(key=lambda path: len(path.name))
        if not candidates:
            raise ValueError(error_msg)
        elif len(candidates) > 1:
            warn(warn_msg)
        return candidates[0]

    def unique_symbols(self, query: ad.AnnData) -> ad.AnnData:
        """Subsetting query to unique feature_name. Renames var_names to feature_name.

        Args:
            query (ad.AnnData): AnnData to subset

        Returns:
            ad.AnnData: Subset AnnData
        """
        query.var["ensembl_id_clean"] = [
            str(val).split(".")[0] for val in query.var["ensembl_id"]
        ]
        subset_query = query.var.copy().dropna(how="any")
        if (
            subset_query[subset_query.duplicated("ensembl_id_clean", keep=False)].shape[
                0
            ]
            != 0
        ):
            print("WARNING: ensembl_id is not unique, need to debug")
        subset_query = subset_query.drop_duplicates("feature_name", keep="first")
        query = query[:, list(subset_query.index)]
        query.var.index = query.var["feature_name"]
        return query

    def into_ref_space(self, query: ad.AnnData, ref: ad.AnnData) -> ad.AnnData:
        """Subsets query to ref.var_names feature space.

        Args:
            query (ad.AnnData): AnnData to subset
            ref (ad.AnnData): Reference feature space

        Returns:
            ad.AnnData: Subset AnnData
        """
        genes_in_common = list(set(ref.var.index).intersection(set(query.var.index)))
        query = query[:, genes_in_common]
        return query

    def get_annotation(self, results: pd.DataFrame) -> pd.DataFrame:
        """Converts frmatch.FRmatch_cell2cluster output into annotation format for adata.obs.
        Applied p-value adjustment and p-value threshold cutoffs.
        Adds "unassigned" for cells with no p-values above cutoff.
        This method is adapted from frmatch.plot_FRmatch_cell2cluster.py.

        Args:
            results (pd.DataFrame): frmatch.FRmatch_cell2cluster output

        Returns:
            pd.DataFrame: Annotation dataframe to be merged with adata.obs.
        """
        p_adj_method = "fdr_by"
        sig_level = 0.05

        # Pivotting results to ref_cluster as index and query_cluster as columns
        results = results.sort_values(
            ["query_cluster", "index", "p_value"], ascending=[True, True, False]
        )
        results = results[
            ~results.duplicated(["index", "query_cluster", "ref_cluster"], keep="first")
        ]
        results = results.pivot(
            index="ref_cluster", columns=["query_cluster", "index"], values="p_value"
        ).replace(np.nan, 0)

        # Adjusting p-values
        results_2 = frmatch.padj_FRmatch(results, p_adj_method=p_adj_method)

        # Setting values to zero unless passes sig_level threshold
        results_3 = results_2.applymap(lambda x: x if x > sig_level else 0)

        # Adding the unassigned row
        samples = results_3.columns[results_3.sum() == 0]
        temp = pd.DataFrame(
            dict(zip(samples, [1] * len(samples))), index=["unassigned"]
        )
        if temp.shape[1] != 0:
            temp.columns.names = ["query_cluster", "index"]
            temp.index.name = "ref_cluster"
            results_4 = pd.concat([results_3, temp])
        else:
            results_4 = results_3.copy()

        # idxmax: returns ref_cluster with largest value per sample
        temp = pd.DataFrame(results_4.max())
        annotation = pd.DataFrame(results_4.idxmax()).reset_index()
        annotation[1] = list(temp[0])

        return annotation

    def copy_annotations(
        self,
        matrix: ad.AnnData,
        annotation: pd.DataFrame,
    ) -> ad.AnnData:
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
        return matrix


def _get_arg_parser():
    parser = add_common_arguments()
    parser.add_argument(
        "--reference-data-dir",
        type=Path,
        required=True,
        help="Path to directory with reference data",
    )
    parser.add_argument(
        "--ensemble-lookup",
        type=Path,
        default="/src/assets/ensemble-lookup.csv",
        help="Ensemble id to gene name csv",
    )
    parser.add_argument("--query-layers-key", help="Data layer to use")
    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    print(parser)
    args = parser.parse_args()
    algorithm = FRmatchAlgorithm()
    result = algorithm.run(**args.__dict__)
    result.save()
