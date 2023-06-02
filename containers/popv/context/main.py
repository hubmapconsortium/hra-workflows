import argparse
from pathlib import Path
from typing import Optional

import numpy as np
import popv
import scanpy
from scanpy import AnnData

TISSUES = [
    "Bladder",
    "Blood",
    "Bone_Marrow",
    "Fat",
    "Heart",
    "Kidney",
    "Large_Intestine",
    "Liver",
    "Lung",
    "Lymph_Node",
    "Mammary",
    "Muscle",
    "Pancreas",
    "Prostate",
    "Salivary Gland",
    "Skin",
    "Small_Intestine",
    "Spleen",
    "Thymus",
    "Trachea",
    "Vasculature",
]


def _get_arg_parser():
    parser = argparse.ArgumentParser(description="Compute annotations using popv")
    parser.add_argument("data", type=scanpy.read_h5ad, help="h5ad data file")
    parser.add_argument("--tissue", choices=TISSUES, required=True, help="Tissue type")
    parser.add_argument(
        "--reference-data",
        type=scanpy.read_h5ad,
        default=AnnData(),
        help="h5ad reference data file",
    )
    parser.add_argument(
        "--models-dir", type=Path, default="./models", help="Models directory"
    )
    parser.add_argument(
        "--cell-ontology-dir",
        type=Path,
        default="./ontology",
        help="Cell ontology directory",
    )
    parser.add_argument("--query-labels-key", help="Data labels key")
    parser.add_argument("--query-batch-key", help="Data batch key")
    parser.add_argument(
        "--ref-labels-key",
        default="cell_ontology_class",
        help="Reference data labels key",
    )
    parser.add_argument(
        "--ref-batch-key", default="donor_assay", help="Reference data batch key"
    )
    parser.add_argument(
        "--unknown-labels-key", default="unknown", help="Unknown labels key"
    )
    parser.add_argument(
        "--samples-per-label", type=int, default=500, help="Number of samples per label"
    )
    parser.add_argument(
        "-n",
        "--sample-size",
        type=int,
        help="Limit data to n randomly selected records",
    )
    parser.add_argument(
        "-o", "--output", type=Path, help="Output file", default="annotations.csv"
    )

    return parser


def _down_sample(data: AnnData, count: int):
    return data[np.random.choice(data.obs_names, count, replace=True)]


def prepare_query(
    data: AnnData,
    reference_data: AnnData,
    *,
    model_dir: Path,
    cell_ontology_dir: Path,
    query_labels_key: Optional[str],
    query_batch_key: Optional[str],
    ref_labels_key: str,
    ref_batch_key: str,
    unknown_labels_key: Optional[str],
    n_samples_per_label: int,
) -> AnnData:
    min_n_samples_per_label = 0
    if ref_labels_key in reference_data.obs.columns:
        min_n_samples_per_label = np.min(
            reference_data.obs.groupby(ref_labels_key).size()
        )

    query = popv.preprocessing.Process_Query(
        data,
        reference_data,
        save_path_trained_models=str(model_dir),
        prediction_mode="fast",
        query_labels_key=query_labels_key,
        query_batch_key=query_batch_key,
        ref_labels_key=ref_labels_key,
        ref_batch_key=ref_batch_key,
        unknown_celltype_label=unknown_labels_key,
        n_samples_per_label=max(n_samples_per_label, min_n_samples_per_label),
        cl_obo_folder=str(cell_ontology_dir) + "/",
        compute_embedding=True,
        hvg=None
        # Maybe set use_gpu=True
    )

    return query.adata


def annotate(query_data: AnnData) -> AnnData:
    popv.annotation.annotate_data(
        query_data,
        # TODO: onclass has been removed due to error in fast mode
        # seen_result_key is not added in fast mode but still expected during compute_consensus
        # https://github.com/YosefLab/PopV/blob/main/popv/annotation.py#L64
        # https://github.com/YosefLab/PopV/blob/main/popv/algorithms/_onclass.py#L199
        methods=["knn_on_scvi", "scanvi", "svm", "rf", "celltypist"]
    )
    return query_data


def main(args: argparse.Namespace):
    data = args.data
    reference_data = args.reference_data
    sample_size = args.sample_size
    if sample_size is not None and sample_size > 0:
        data = _down_sample(data, sample_size)
        reference_data = _down_sample(reference_data, sample_size)

    query_data = prepare_query(
        data,
        reference_data,
        model_dir=args.models_dir / f"pretrained_models_{args.tissue}_ts",
        cell_ontology_dir=args.cell_ontology_dir,
        query_labels_key=args.query_labels_key,
        query_batch_key=args.query_batch_key,
        ref_labels_key=args.ref_labels_key,
        ref_batch_key=args.ref_batch_key,
        unknown_labels_key=args.unknown_labels_key,
        n_samples_per_label=args.samples_per_label,
    )

    result = annotate(query_data)
    result.obs.to_csv(args.output)


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)