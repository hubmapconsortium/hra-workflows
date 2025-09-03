import csv
import logging
import typing as t
from pathlib import Path

from panhumanpy import AzimuthNN
import pandas
import scanpy

from src.algorithm import Algorithm, RunResult, add_common_arguments
from src.util.layers import set_data_layer


class PanHumanAzimuthOrganMetadata(t.TypedDict):
    model: str


class PanHumanAzimuthOptions(t.TypedDict):
    ensemble_lookup: Path
    query_layers_key: t.Optional[str]


class PanHumanAzimuthAlgorithm(Algorithm[PanHumanAzimuthOrganMetadata, PanHumanAzimuthOptions]):
    def __init__(self):
        super().__init__("azimuth_fine", is_pan_organ=True)

    def do_run(
        self,
        matrix: Path,
        organ: str,
        metadata: PanHumanAzimuthOrganMetadata,
        options: PanHumanAzimuthOptions,
    ) -> RunResult:
        """Annotate data using pan-human-azimuth."""
        data = scanpy.read_h5ad(matrix)
        data = set_data_layer(data, options["query_layers_key"])
        data = self.clean(data)
        data = self.normalize(data)
        data, var_names = self.normalize_var_names(data, options)

        azimuth = AzimuthNN(data,model_version="v1")
        azimuth.azimuth_refine()
        data = azimuth.pack_adata()

        data.var_names = t.cast(t.Any, var_names)

        return {"data": data, "organ_id": "UBERON:0013702", "organ_level": "body_proper"}

    def clean(self, data: scanpy.AnnData) -> scanpy.AnnData:
        """Cleans the data removing any incompatible preprocessing that may exist.

        Args:
            data (scanpy.AnnData): Original data to clean

        Returns:
            scanpy.AnnData: Clean data
        """
        data.obsm = None
        return data

    def normalize(self, data: scanpy.AnnData) -> scanpy.AnnData:
        """Normalizes data according to pan-human-azimuth requirements.

        Args:
            data (scanpy.AnnData): Original data to be normalized

        Returns:
            scanpy.AnnData: Normalized data
        """
        primary_column = "feature_name"
        alternative_primary_columns = ["hugo_symbol", "gene_symbol"]
        if primary_column not in data.var.columns:
            logging.warning(
                f"Missing {primary_column} data column. Attempting to use {alternative_primary_columns} instead"
            )
            for column in alternative_primary_columns:
                if column in data.var.columns:
                    data.var = data.var.rename(columns={column: primary_column})
                    break

        data.var_names_make_unique()
        return data

    def normalize_var_names(
        self, data: scanpy.AnnData, options: PanHumanAzimuthOptions
    ) -> t.Tuple[scanpy.AnnData, pandas.Index]:
        """Normalizes variable names, replacing ensemble ids with the corresponding gene name.

        Args:
            data (scanpy.AnnData): Data with potentially non-normalized names
            options (PanHumanAzimuthOptions): Options containing the ensemble id mapping file path

        Returns:
            t.Tuple[scanpy.AnnData, pandas.Index]: The normalized data along with the original names
        """
        lookup = self.load_ensemble_lookup(options)
        names = data.var_names

        def getNewName(name: str):
            key = name.split(".", 1)[0]
            return lookup.get(key, name)

        data.var_names = t.cast(t.Any, names.map(getNewName))
        return data, names

    def load_ensemble_lookup(self, options: PanHumanAzimuthOptions) -> t.Dict[str, str]:
        """Load a file mapping ensemble id to gene names.

        Args:
            options (PanHumanAzimuthOptions): Options with the mapping file path

        Returns:
            t.Dict[str, str]: Loaded mapping
        """
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
        default="/src/assets/ensemble-lookup.csv",
        help="Ensemble id to gene name csv",
    )
    parser.add_argument("--query-layers-key", help="Data layer to use")

    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    algorithm = PanHumanAzimuthAlgorithm()
    result = algorithm.run(**args.__dict__)
    result.save()
