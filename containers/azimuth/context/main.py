import logging
import subprocess
import typing as t
from pathlib import Path

import anndata
import pandas

from src.algorithm import Algorithm, OrganLookup, add_common_arguments


class AzimuthOptions(t.TypedDict):
    reference_data_dir: Path


class AzimuthOrganLookup(OrganLookup[str]):
    def __init__(self, mapping_file: Path):
        super().__init__(mapping_file)

    def get_builtin_options(self):
        # TODO read from mapping file?
        references = [
            "adiposeref",
            "bonemarrowref",
            "fetusref",
            "heartref",
            "humancortexref",
            "kidneyref",
            "lungref",
            "mousecortexref",
            "pancreasref",
            "pbmcref",
            "tonsilref",
        ]
        return zip(references, references)


class AzimuthAlgorithm(Algorithm[str, AzimuthOptions]):
    def __init__(self):
        super().__init__(AzimuthOrganLookup, "predicted.ann_finest_level")

    def do_run(self, matrix: Path, organ: str, options: AzimuthOptions):
        data = anndata.read_h5ad(matrix)
        reference_data = self.find_reference_data(organ, options["reference_data_dir"])

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

        return data

    def create_clean_matrix(self, matrix: anndata.AnnData):
        clean_obs = pandas.DataFrame(index=matrix.obs.index)
        clean_matrix = matrix.copy()
        clean_matrix.obs = clean_obs
        return clean_matrix

    def copy_annotations(
        self, matrix: anndata.AnnData, annotated_matrix: anndata.AnnData
    ):
        matrix.obs = matrix.obs.join(annotated_matrix.obs, rsuffix="_azimuth")

    def run_azimuth_scripts(self, matrix_path: Path, reference_data: Path):
        script_command = ["Rscript", "/run_azimuth.R", matrix_path, reference_data]
        subprocess.run(script_command, capture_output=True, check=True, text=True)
        return "./result.h5ad"

    def find_reference_data(self, organ: str, dir: Path):
        def is_reference_data_candidate(path: Path):
            return path.is_dir() and organ.lower() in path.name.lower()

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
    ):
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
