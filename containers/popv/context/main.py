import typing as t
from logging import warn
from pathlib import Path

import numpy
import popv
import scanpy
import anndata
import torch
import pandas

# import typing as t
import csv

from src.algorithm import Algorithm, OrganLookup, add_common_arguments


class PopvOptions(t.TypedDict):
    reference_data_dir: Path
    models_dir: Path
    prediction_mode: str
    cell_ontology_dir: str
    query_labels_key: t.Optional[str]
    query_batch_key: t.Optional[str]
    query_layers_key: t.Optional[str]
    ref_labels_key: str
    ref_batch_key: str
    unknown_labels_key: str
    samples_per_label: int
    data_genes_column: str
    ensemble_lookup: Path


class PopvAlgorithm(Algorithm[str, PopvOptions]):
    def __init__(self):
        super().__init__(OrganLookup)
        # self.models_dir = Path(os.environ["MODELS_DIR"])

    def do_run(self, matrix: Path, organ: str, options: PopvOptions):
        data = scanpy.read_h5ad(matrix)
        data = self.prepare_query(data, organ, options)
        popv.annotation.annotate_data(
            data,
            # TODO: onclass has been removed due to error in fast mode
            # seen_result_key is not added to the result in fast mode but still expected during compute_consensus
            # https://github.com/YosefLab/PopV/blob/main/popv/annotation.py#L64
            # https://github.com/YosefLab/PopV/blob/main/popv/algorithms/_onclass.py#L199
            methods=["knn_on_scvi", "scanvi", "svm", "rf", "celltypist"],
        )
        return data

    def prepare_query(
        self, data: scanpy.AnnData, organ: str, options: PopvOptions
    ) -> scanpy.AnnData:
        reference_data_path = self.find_reference_data(
            options["reference_data_dir"], organ
        )
        model_path = self.find_model_dir(options["models_dir"], organ)
        reference_data = scanpy.read_h5ad(reference_data_path)
        data, var_names = self.normalize_var_names(data, options)
        # data.var_names = t.cast(t.Any, var_names)
        updated_data = self.add_model_genes(
            data,
            model_path,
            options["data_genes_column"],
        )
        n_samples_per_label = self.get_n_samples_per_label(reference_data, options)
        query = popv.preprocessing.Process_Query(
            updated_data,
            reference_data,
            save_path_trained_models=str(model_path),
            prediction_mode=options["prediction_mode"],
            query_labels_key=options["query_labels_key"],
            query_batch_key=options["query_batch_key"],
            query_layers_key=options["query_layers_key"],
            ref_labels_key=options["ref_labels_key"],
            ref_batch_key=options["ref_batch_key"],
            unknown_celltype_label=options["unknown_labels_key"],
            n_samples_per_label=n_samples_per_label,
            cl_obo_folder=f"{options['cell_ontology_dir']}/",
            compute_embedding=True,
            hvg=None,
            use_gpu=False,  # Using gpu with docker requires additional setup
        )

        return query.adata

    def get_n_samples_per_label(
        self, reference_data: scanpy.AnnData, options: PopvOptions
    ) -> int:
        ref_labels_key = options["ref_labels_key"]
        n_samples_per_label = options["samples_per_label"]
        if ref_labels_key in reference_data.obs.columns:
            n = numpy.min(reference_data.obs.groupby(ref_labels_key).size())
            n_samples_per_label = numpy.max((n_samples_per_label, t.cast(int, n)))
        return n_samples_per_label

    def find_reference_data(self, dir: Path, organ: str) -> Path:
        def is_reference_data_candidate(path: Path):
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

    def find_model_dir(self, dir: Path, organ: str) -> Path:
        def is_model_candidate(path: Path):
            return path.is_dir() and organ.lower() in path.name.lower()

        return self._find_in_dir(
            dir,
            is_model_candidate,
            f"Cannot find model directory for organ '{organ}'",
            f"Multiple model directory candidates for organ '{organ}'",
        )

    def _find_in_dir(
        self, dir: Path, cond: t.Callable[[Path], bool], error_msg: str, warn_msg: str
    ):
        candidates = list(filter(cond, dir.iterdir()))
        candidates.sort(key=lambda path: len(path.name))

        if not candidates:
            raise ValueError(error_msg)
        elif len(candidates) > 1:
            warn(warn_msg)
        return candidates[0]

    def normalize_var_names(
        self, data: scanpy.AnnData, options: PopvOptions
    ) -> t.Tuple[scanpy.AnnData, pandas.Index]:
        lookup = self.load_ensemble_lookup(options)
        names = data.var_names

        def getNewName(name: str):
            key = name.split(".", 1)[0]
            return lookup.get(key, name)

        data.var_names = t.cast(t.Any, names.map(getNewName))
        return data, names

    def load_ensemble_lookup(self, options: PopvOptions):
        with open(options["ensemble_lookup"]) as file:
            reader = csv.DictReader(file)
            lookup: t.Dict[str, str] = {}
            for row in reader:
                lookup[row["ensemble"]] = row["gene_name"]
        return lookup

    def add_model_genes(
        self,
        data: scanpy.AnnData,
        model_path: Path,
        data_genes_column: str,
    ) -> scanpy.AnnData:
        """Adds genes from model not present in input data to input data. Needed for preprocessing bug"""
        model_genes = torch.load(
            Path.joinpath(model_path, "scvi/model.pt"), map_location="cpu"
        )["var_names"]
        n_obs_data = data.X.shape[0]
        new_genes = set(numpy.setdiff1d(model_genes, data.var_names))
        new_data = scanpy.AnnData(
            X=numpy.zeros((n_obs_data, len(new_genes))), var=new_genes
        )
        new_data.obs_names = data.obs_names
        new_data.var_names = new_genes
        merged_data = anndata.concat([data, new_data], axis=1)
        print(merged_data.shape)
        merged_data_clean = self.combine_duplicate_columns(merged_data)
        print(merged_data_clean.shape)
        return merged_data

    def combine_duplicate_columns(self, data: scanpy.AnnData) -> scanpy.AnnData:
        """Combines columns with same gene name takes their mean"""
        # INCOMPLETE
        duplicated_var_names = data.var_names[
            data.var_names.duplicated(keep="first")
        ].unique()
        averaged_data = []
        for duplicated_var_name in duplicated_var_names:
            averaged_data.append(data[:, duplicated_var_name].X.mean(axis=1))
        new_data = anndata.AnnData(
            X=averaged_data, obs=data.obs, var=duplicated_var_names
        )
        data = anndata.concat(
            [data[:, ~data.var_names.isin(duplicated_var_names)], new_data]
        )
        return data


def _get_arg_parser():
    parser = add_common_arguments()
    parser.add_argument(
        "--reference-data-dir",
        type=Path,
        required=True,
        help="Path to directory with reference data",
    )
    parser.add_argument(
        "--models-dir",
        type=Path,
        required=True,
        help="Path to models directory",
    )
    parser.add_argument(
        "--data-genes-column", required=True, help="Data gene names column"
    )
    parser.add_argument("--prediction-mode", default="fast", help="Prediction mode")
    parser.add_argument(
        "--cell-ontology-dir",
        type=Path,
        default="/ontology",
        help="Cell ontology directory",
    )
    parser.add_argument("--query-labels-key", help="Data labels key")
    parser.add_argument("--query-batch-key", help="Data batch key")
    parser.add_argument("--query-layers-key", help="Name of layer with raw counts")
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
        "--ensemble-lookup",
        type=Path,
        default="/ensemble-lookup.csv",
        help="Ensemble id to gene name csv",
    )
    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    algorithm = PopvAlgorithm()
    result = algorithm.run(**args.__dict__)
    result.save()
