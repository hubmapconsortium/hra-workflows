import subprocess
from pathlib import Path

import anndata
import pandas

from src.algorithm import Algorithm, OrganLookup, add_common_arguments


class AzimuthOrganLookup(OrganLookup[str]):
    def __init__(self, mapping_file: Path):
        super().__init__(mapping_file)

    def get_builtin_options(self):
        references = ["adiposeref" ,"bonemarrowref" ,"fetusref" ,"heartref"  ,"humancortexref"   ,"kidneyref"  ,"lungref"   ,"mousecortexref"   ,"pancreasref"   ,"pbmcref"  ,"tonsilref"]
        return zip(references, references)


class AzimuthAlgorithm(Algorithm[str, dict]):
    def __init__(self):
        super().__init__(AzimuthOrganLookup)

    def do_run(self, matrix: Path, organ: str, options: dict):
        print('foo')
        data = anndata.read_h5ad(matrix)

        # Azimuth chokes when trying to load matrices that has
        # obs columns of dtype 'object'. As a workaround we create a
        # clean matrix without obs columns on which azimuth is run
        # after which the annotations are copied back to the original matrix
        clean_matrix_path = Path("clean_matrix.h5ad")
        clean_matrix = self.create_clean_matrix(data)
        clean_matrix.write_h5ad(clean_matrix_path)
        annotated_matrix_path = self.run_azimuth_scripts(clean_matrix_path, organ)
        annotated_matrix = anndata.read_h5ad(annotated_matrix_path)
        self.copy_annotations(data, annotated_matrix)

        return annotated_matrix

    def create_clean_matrix(self, matrix: anndata.AnnData):
        clean_obs = pandas.DataFrame(index=matrix.obs.index)
        clean_matrix = matrix.copy()
        clean_matrix.obs = clean_obs
        return clean_matrix

    def copy_annotations(
        self, matrix: anndata.AnnData, annotated_matrix: anndata.AnnData
    ):
        matrix.obs = matrix.obs.join(annotated_matrix.obs)

    def run_azimuth_scripts(self, matrix_path: Path, organ: str):
        script_outputs = [
            "./result.h5ad"
            # "./version_metadata.json",
            # "./annotations.csv",
        ]
        script_command = [
            "Rscript",
            "run_azimuth.R",
            organ,
            matrix_path
        ]
        subprocess.run(script_command, capture_output=False, check=True)
        return script_outputs[0]


if __name__ == "__main__":
    parser = add_common_arguments()
    args = parser.parse_args()
    algorithm = AzimuthAlgorithm()
    result = algorithm.run(**args.__dict__)
    result.save()
