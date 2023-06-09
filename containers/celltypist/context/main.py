import argparse
import itertools
import json
import logging
from pathlib import Path

import celltypist
import scanpy
from scanpy import AnnData


def _organ_or_model(value: str):
    with open("/organ-mapping.json") as mapping_file:
        mapping = json.load(mapping_file)

    value = value.lower()
    models = celltypist.models.get_all_models()
    items = itertools.chain(mapping.values(), zip(models, models))
    for key, model in items:
        if key.lower() == value:
            return celltypist.models.Model.load(model)

    raise ValueError("Invalid organ")


def _get_arg_parser():
    parser = argparse.ArgumentParser(description="Compute annotations using celltypist")
    parser.add_argument("data", type=scanpy.read_h5ad, help="h5ad data file")
    parser.add_argument(
        "--organ",
        type=_organ_or_model,
        required=True,
        dest="model",
        help="Organ uberon id",
    )
    parser.add_argument(
        "--existing-annotations-column", help="Column with existing annotations"
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


def normalize(data: AnnData) -> AnnData:
    primary_column = "feature_name"
    alternative_primary_column = "gene_symbol"
    if primary_column not in data.var.columns:
        logging.warning(
            f"Missing {primary_column} data column. Attempting to use {alternative_primary_column} instead"
        )
        data.var = data.var.rename({alternative_primary_column: primary_column})

    # Note: I've excluded some of the original code from an if statement that is always false
    # https://github.com/cns-iu/ct-ann-predictive-analytics/blob/main/celltypist/celltypist_prediction_pipeline.py#L116

    scanpy.pp.normalize_total(data, target_sum=1e4)
    scanpy.pp.log1p(data)
    data.var_names_make_unique()
    return data


def annotate(data: AnnData, model: celltypist.models.Model) -> AnnData:
    return celltypist.annotate(data, model, majority_voting=True).to_adata()


def main(args: argparse.Namespace):
    data = normalize(args.data)
    annotations = annotate(data, args.model)
    existing_annotations_column = args.existing_annotations_column

    if existing_annotations_column:
        matches = lambda row: row[existing_annotations_column] == row["majority_voting"]
        annotations.obs["exp_vs_pred"] = annotations.obs.apply(matches, axis=1)

    annotations.write_h5ad(args.output_matrix)
    annotations.obs.to_csv(args.output, index=True)


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
