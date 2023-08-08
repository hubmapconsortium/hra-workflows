import logging
from pathlib import Path

import celltypist
import scanpy

# From shared docker image
from shared.algorithm import Algorithm, OrganLookup, add_common_arguments


class CelltypistOrganLookup(OrganLookup[celltypist.Model]):
    def get_builtin_options(self):
        models = celltypist.models.get_all_models()
        return map(lambda model: (model, self.from_raw(model)), models)

    def from_raw(self, id: str):
        return celltypist.models.Model.load(id)


class CelltypistAlgorithm(Algorithm[celltypist.Model, dict]):
    def __init__(self):
        super().__init__(CelltypistOrganLookup)

    def do_run(self, matrix: Path, organ: celltypist.Model, options: dict):
        data = scanpy.read_h5ad(matrix)
        data = self.normalize(data)
        return celltypist.annotate(data, organ, majority_voting=True).to_adata()

    def normalize(self, data: scanpy.AnnData) -> scanpy.AnnData:
        primary_column = "feature_name"
        alternative_primary_column = "gene_symbol"
        if primary_column not in data.var.columns:
            logging.warning(
                f"Missing {primary_column} data column. Attempting to use {alternative_primary_column} instead"
            )
            data.var = data.var.rename({alternative_primary_column: primary_column})

        # NOTE: I've excluded some of the original code from an if statement that is always false
        # https://github.com/cns-iu/ct-ann-predictive-analytics/blob/main/celltypist/celltypist_prediction_pipeline.py#L116

        scanpy.pp.normalize_total(data, target_sum=1e4)
        scanpy.pp.log1p(data)
        data.var_names_make_unique()
        return data


if __name__ == "__main__":
    parser = add_common_arguments()
    args = parser.parse_args()
    algorithm = CelltypistAlgorithm()
    result = algorithm.run(**args.__dict__)
    result.save()
