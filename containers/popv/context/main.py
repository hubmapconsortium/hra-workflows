import os
import typing as t
from pathlib import Path

import numpy as np
import popv
import scanpy

# From shared docker image
from shared.algorithm import Algorithm, OrganLookup, add_common_arguments


class PopvOptions(t.TypedDict):
    reference_data: Path
    prediction_mode: str
    cell_ontology_dir: str
    query_labels_key: t.Optional[str]
    query_batch_key: t.Optional[str]
    query_layers_key: t.Optional[str]
    ref_labels_key: str
    ref_batch_key: str
    unknown_labels_key: str
    samples_per_label: int


class PopvOrganLookup(OrganLookup[str]):
    def __init__(self, mapping_file: Path):
        super().__init__(mapping_file)
        self.models_dir = Path(os.environ["MODELS_DIR"])
        self.prefix = "pretrained_models_"
        self.suffix = "_ts"

    def get_builtin_options(self):
        dirs = self.get_dirs()
        return map(lambda dir: (self.get_name(dir), dir), dirs)

    def get_dirs(self) -> t.Iterable[Path]:
        return self.models_dir.glob(f"{self.prefix}*{self.suffix}")

    def get_name(self, dir: Path) -> str:
        start = len(self.prefix)
        end = -len(self.suffix)
        return dir.name[start:end]


class PopvAlgorithm(Algorithm[str, PopvOptions]):
    def __init__(self):
        super().__init__(PopvOrganLookup)
        self.models_dir = Path(os.environ["MODELS_DIR"])

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
        query = popv.preprocessing.Process_Query(
            data,
            options["reference_data"],
            save_path_trained_models=str(self.models_dir / organ),
            prediction_mode=options["prediction_mode"],
            query_labels_key=options["query_labels_key"],
            query_batch_key=options["query_batch_key"],
            query_layers_key=options["query_layers_key"],
            ref_labels_key=options["ref_labels_key"],
            ref_batch_key=options["ref_batch_key"],
            unknown_celltype_label=options["unknown_labels_key"],
            n_samples_per_label=self.get_n_samples_per_label(options),
            cl_obo_folder=f"{options['cell_ontology_dir']}/",
            compute_embedding=True,
            hvg=None,
            use_gpu=False,  # Using gpu with docker requires additional setup
        )

        return query.adata

    def get_n_samples_per_label(self, options: PopvOptions) -> int:
        reference_data = options["reference_data"]
        ref_labels_key = options["ref_labels_key"]
        n_samples_per_label = options["samples_per_label"]
        if ref_labels_key in reference_data.obs.columns:
            n = np.min(reference_data.obs.groupby(ref_labels_key).size())
            n_samples_per_label = max(n_samples_per_label, n)
        return n_samples_per_label


def _get_arg_parser():
    parser = add_common_arguments()
    parser.add_argument(
        "--reference-data",
        type=scanpy.read_h5ad,
        required=True,
        help="h5ad reference data file",
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

    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    algorithm = PopvAlgorithm()
    result = algorithm.run(**args.__dict__)
    result.save()
