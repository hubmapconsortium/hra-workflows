import logging
import subprocess
import typing as t
from pathlib import Path

import anndata
import pandas

from src.algorithm import Algorithm, RunResult, add_common_arguments


class AzimuthOrganMetadata(t.TypedDict):
    model: str
    organ_level: str
    prediction_column: str


class AzimuthOptions(t.TypedDict):
    reference_data_dir: Path


class AzimuthAlgorithm(Algorithm[AzimuthOrganMetadata, AzimuthOptions]):
    def do_run(
        self,
        matrix: Path,
        organ: str,
        metadata: AzimuthOrganMetadata,
        options: AzimuthOptions,
    ) -> RunResult:
        """Annotate data using azimuth."""
        data = anndata.read_h5ad(matrix)
        reference_data = self.find_reference_data(
            organ, metadata["model"], options["reference_data_dir"]
        )

        # Azimuth chokes when trying to load matrices that has
        # obs columns of dtype 'object'. As a workaround we create a
        # clean matrix without obs columns on which azimuth is run
        # after which the annotations are copied back to the original matrix
        clean_matrix_path = Path("clean_matrix.h5ad")
        clean_matrix = self.create_clean_matrix(data)
        clean_matrix.write_h5ad(clean_matrix_path)

        annotated_matrix_path = self.run_azimuth_scripts(
            clean_matrix_path, reference_data
        )
        annotated_matrix = anndata.read_h5ad(annotated_matrix_path)
        self.copy_annotations(data, annotated_matrix)

        return {
            "data": data,
            "organ_level": metadata["organ_level"],
            "prediction_column": "predicted." + metadata["prediction_column"],
        }

    def create_clean_matrix(self, matrix: anndata.AnnData) -> anndata.AnnData:
        """Creates a copy of the data with all observation columns removed.

        Args:
            matrix (anndata.AnnData): Original data

        Returns:
            anndata.AnnData: Cleaned data
        """
        clean_obs = pandas.DataFrame(index=matrix.obs.index)
        clean_matrix = matrix.copy()
        clean_matrix.obs = clean_obs
        return clean_matrix

    def copy_annotations(
        self, matrix: anndata.AnnData, annotated_matrix: anndata.AnnData
    ) -> None:
        """Copies annotations from one matrix to another.

        Args:
            matrix (anndata.AnnData): Matrix to copy to
            annotated_matrix (anndata.AnnData): Matrix to copy from
        """
        matrix.obs = matrix.obs.join(annotated_matrix.obs, rsuffix="_azimuth")

    def run_azimuth_scripts(self, matrix_path: Path, reference_data: Path) -> str:
        """Creates a subprocess running the Azimuth annotation R script.

        Args:
            matrix_path (Path): Path to data file
            reference_data (Path): Path to model reference data directory

        Returns:
            str: Path to the output data file
        """
        script_command = ["Rscript", "/run_azimuth.R", matrix_path, reference_data]
        subprocess.run(script_command, capture_output=True, check=True, text=True)
        return "./result.h5ad"

    def find_reference_data(self, organ: str, model: str, dir: Path) -> Path:
        """Finds the reference data directory for a model.

        Args:
            organ (str): Organ id
            model (str): Model name
            dir (Path): Directory to search

        Raises:
            ValueError: If no reference data could be found

        Returns:
            Path: The data directory
        """

        def is_reference_data_candidate(path: Path):
            return path.is_dir() and model.lower() in path.name.lower()

        subdir = self._find_in_dir(
            dir,
            is_reference_data_candidate,
            f"Cannot find reference data for organ '{organ}'",
            f"Multiple reference data candidates for organ '{organ}'",
        )
        # idx.annoy and ref.Rds is always located inside an 'azimuth' subdirectory
        return subdir / "azimuth"

    def _find_in_dir(
        self, dir: Path, cond: t.Callable[[Path], bool], error_msg: str, warn_msg: str
    ) -> Path:
        """Search a directory for a entry which passes the provided test.

        Args:
            dir (Path): Directory to search
            cond (t.Callable[[Path], bool]): Test used to match sub entries
            error_msg (str): Error message used when no entries match
            warn_msg (str): Warning message use when multiple entries match

        Raises:
            ValueError: If there are no matching sub entries

        Returns:
            Path:
                The matching entry.
                If multiple entries match the one with the shortest name is returned.
        """
        candidates = list(filter(cond, dir.iterdir()))
        candidates.sort(key=lambda path: len(path.name))

        if not candidates:
            raise ValueError(error_msg)
        elif len(candidates) > 1:
            logging.warn(warn_msg)
        return candidates[0]


def _get_arg_parser():
    parser = add_common_arguments()
    parser.add_argument(
        "--reference-data-dir",
        type=Path,
        required=True,
        help="Path to directory with reference data",
    )

    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    algorithm = AzimuthAlgorithm()
    result = algorithm.run(**args.__dict__)
    result.save()
