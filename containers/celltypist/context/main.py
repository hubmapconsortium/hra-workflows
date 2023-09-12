import csv
import logging
import typing as t
from pathlib import Path

import celltypist
import pandas
import scanpy

from src.algorithm import Algorithm, OrganLookup, add_common_arguments


class CelltypistOptions(t.TypedDict):
    ensemble_lookup: Path


class CelltypistOrganLookup(OrganLookup[celltypist.Model]):
    def __init__(self, mapping_file: Path):
        super().__init__(mapping_file)

    def get_builtin_options(self):
        models = celltypist.models.get_all_models()
        return map(lambda model: (model, self.from_raw(model)), models)

    def from_raw(self, id: str):
        return celltypist.models.Model.load(id)


class CelltypistAlgorithm(Algorithm[celltypist.Model, CelltypistOptions]):
    def __init__(self):
        super().__init__(CelltypistOrganLookup)

    def do_run(self, matrix: Path, organ: celltypist.Model, options: CelltypistOptions):
        data = scanpy.read_h5ad(matrix)
        data = self.normalize(data)
        data, var_names = self.normalize_var_names(data, options)
        data = celltypist.annotate(data, organ, majority_voting=True).to_adata()
        data.var_names = t.cast(t.Any, var_names)
        return data

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

    def normalize_var_names(
        self, data: scanpy.AnnData, options: CelltypistOptions
    ) -> t.Tuple[scanpy.AnnData, pandas.Index]:
        lookup = self.load_ensemble_lookup(options)
        names = data.var_names

        def getNewName(name: str):
            key = name.split(".", 1)[0]
            return lookup.get(key, name)

        data.var_names = t.cast(t.Any, names.map(getNewName))
        return data, names

    def load_ensemble_lookup(self, options: CelltypistOptions):
        with open(options["ensemble_lookup"]) as file:
            reader = csv.DictReader(file)
            lookup: t.Dict[str, str] = {}
            for row in reader:
                lookup[row["ensemble"]] = row["gene_name"]
        return lookup


def _get_arg_parser():
    parser = add_common_arguments()
    parser.add_argument(
        "--ensemble-lookup",
        type=Path,
        default="/ensemble-lookup.csv",
        help="Ensemble id to gene name csv",
    )

    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    algorithm = CelltypistAlgorithm()
    result = algorithm.run(**args.__dict__)
    result.save()
