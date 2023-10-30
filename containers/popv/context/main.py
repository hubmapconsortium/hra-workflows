import typing as t
from logging import warn
from pathlib import Path

import numpy
import popv
import scanpy

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
    ref_gene_column: str
    data_gene_column: str


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
        filtered_data = self.filter_genes(
            data,
            reference_data,
            options["ref_gene_column"],
            options["data_gene_column"],
        )
        n_samples_per_label = self.get_n_samples_per_label(reference_data, options)

        query = popv.preprocessing.Process_Query(
            filtered_data,
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

    def filter_genes(
        data: scanpy.AnnData,
        reference_data: scanpy.AnnData,
        ref_gene_column: str,
        data_gene_column: str,
    ) -> scanpy.AnnData:
        """PoPV preprocessing fails on encountering non reference data genes. Filters data to have only reference data genes"""
        reference_data_genes = reference_data.var[ref_gene_column].unique().tolist()
        filtered_data_var = data.var[
            data.var[data_gene_column].isin(reference_data_genes)
        ]
        numerical_data_var_index = data.var.index.get_indexer(
            filtered_data_var.index.tolist()
        )
        return data[:, numerical_data_var_index]


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
        "--ref-genes-column",
        default="feature_name",
        help="Reference data gene names column",
    )
    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    algorithm = PopvAlgorithm()
    result = algorithm.run(**args.__dict__)
    result.save()
