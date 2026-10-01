import logging
import typing as t
from pathlib import Path

import scanpy

from src.algorithm import Algorithm, RunResult, add_common_arguments
from src.util.layers import set_data_layer

DEFAULT_ORGAN_LEVEL = "author"


class AuthorOrganMetadata(t.TypedDict):
    pass


class AuthorOptions(t.TypedDict):
    label_column: str
    organ_level: t.Optional[str]
    query_layers_key: t.Optional[str]
    keep_unlabeled: bool


class AuthorAlgorithm(Algorithm[AuthorOrganMetadata, AuthorOptions]):
    """Uses the cell types assigned by the dataset authors as annotations."""

    def __init__(self):
        super().__init__(is_pan_organ=True)

    def do_run(
        self,
        matrix: Path,
        organ: str,
        metadata: AuthorOrganMetadata,
        options: AuthorOptions,
    ) -> RunResult:
        label_column = options["label_column"]
        data = scanpy.read_h5ad(matrix)
        if label_column not in data.obs.columns:
            raise ValueError(f"Missing author cell type column '{label_column}'")

        if not options["keep_unlabeled"]:
            data = self.remove_unlabeled(data, label_column)
        if data.n_obs == 0:
            raise ValueError("No cells with an author assigned cell type")

        data = set_data_layer(data, options["query_layers_key"])
        data = self.normalize(data)
        data.obs[label_column] = data.obs[label_column].astype(str)

        return {
            "data": data,
            "organ_level": options["organ_level"] or DEFAULT_ORGAN_LEVEL,
            "prediction_column": label_column,
        }

    def remove_unlabeled(self, data: scanpy.AnnData, label_column: str) -> scanpy.AnnData:
        """Removes cells without an author assigned cell type.

        Args:
            data (scanpy.AnnData): Original data
            label_column (str): Column with the author labels

        Returns:
            scanpy.AnnData: Data containing only labeled cells
        """
        labels = data.obs[label_column]
        mask = labels.notna() & ~labels.astype(str).str.strip().str.lower().isin(["", "nan", "none"])
        logging.info(f"Removing {(~mask).sum()} cells without an author assigned cell type")
        return data[mask].copy()

    def normalize(self, data: scanpy.AnnData) -> scanpy.AnnData:
        """Log1p normalizes data with 10,000 counts per cell for gene expression calculations.

        Args:
            data (scanpy.AnnData): Original data to be normalized

        Returns:
            scanpy.AnnData: Normalized data
        """
        data.obsm = None
        scanpy.pp.normalize_total(data, target_sum=1e4)
        scanpy.pp.log1p(data)
        data.var_names_make_unique()
        return data


def _get_arg_parser():
    parser = add_common_arguments()
    parser.add_argument("--label-column", required=True, help="Obs column with author cell types")
    parser.add_argument("--organ-level", help="Organ level used for crosswalking")
    parser.add_argument("--query-layers-key", help="Data layer to use")
    parser.add_argument(
        "--keep-unlabeled",
        action="store_true",
        help="Keep cells without an author assigned cell type",
    )

    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    algorithm = AuthorAlgorithm()
    result = algorithm.run(**args.__dict__)
    result.save()
