import argparse
import itertools
import json
import os
from pathlib import Path
from typing import Optional

import numpy as np
import popv
import scanpy
from scanpy import AnnData

_MODELS_DIR = Path(os.environ["MODELS_DIR"])


def _down_sample(data: AnnData, count: int):
    return data[np.random.choice(data.obs_names, count, replace=True)]


def _find_models():
    prefix = "pretrained_models_"
    suffix = "_ts"
    dirs = _MODELS_DIR.glob(f"{prefix}*{suffix}")
    return [(d.name[len(prefix) : -len(suffix)], d.name) for d in dirs]


def _organ_or_tissue(value: str):
    with open("/organ-mapping.json") as mapping_file:
        mapping = json.load(mapping_file)

    value = value.lower()
    items = itertools.chain(mapping.values(), _find_models())
    for key, tissue in items:
        if key.lower() == value:
            return tissue

    raise ValueError("Invalid organ")


def _get_arg_parser():
    parser = argparse.ArgumentParser(description="Compute annotations using popv")
    parser.add_argument("data", type=scanpy.read_h5ad, help="h5ad data file")
    parser.add_argument(
        "--organ",
        type=_organ_or_tissue,
        required=True,
        dest="tissue",
        help="Organ uberon id",
    )
    parser.add_argument(
        "--reference-data",
        type=scanpy.read_h5ad,
        default=AnnData(),
        help="h5ad reference data file",
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
        "-o", "--output", type=Path, default="annotations.csv", help="Output file"
    )
    parser.add_argument(
        "--output-matrix",
        type=Path,
        default="annotated_matrix.h5ad",
        help="Annotated matrix output file",
    )

    return parser


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
        methods=["knn_on_scvi", "scanvi", "svm", "rf", "celltypist"],
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
        model_dir=_MODELS_DIR / args.tissue,
        cell_ontology_dir=args.cell_ontology_dir,
        query_labels_key=args.query_labels_key,
        query_batch_key=args.query_batch_key,
        ref_labels_key=args.ref_labels_key,
        ref_batch_key=args.ref_batch_key,
        unknown_labels_key=args.unknown_labels_key,
        n_samples_per_label=args.samples_per_label,
    )

    result = annotate(query_data)
    result.write_h5ad(args.output_matrix)
    result.obs.to_csv(args.output)


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
