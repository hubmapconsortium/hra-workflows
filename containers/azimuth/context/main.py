import subprocess
from pathlib import Path

import write_metadata  # type: ignore From azimuth-annotate docker image

from src.algorithm import Algorithm, OrganLookup, add_common_arguments


class AzimuthOrganLookup(OrganLookup[str]):
    def __init__(self, mapping_file: Path):
        super().__init__(mapping_file)

    def get_builtin_options(self):
        references = ["RK", "LK", "RL", "LL", "HT"]
        return zip(references, references)


class AzimuthAlgorithm(Algorithm[str, dict]):
    def __init__(self):
        super().__init__(AzimuthOrganLookup)

    def do_run(self, matrix: Path, organ: str, options: dict):
        script_outputs = [
            "./secondary_analysis.h5ad",
            "./version_metadata.json",
            "./annotations.csv",
        ]
        script_command = [
            "Rscript",
            "/azimuth_analysis.R",
            organ,
            matrix,
            matrix,
        ]

        subprocess.run(script_command, capture_output=True, check=True)
        write_metadata.main(*script_outputs)
        return Path(script_outputs[0])


if __name__ == "__main__":
    parser = add_common_arguments()
    args = parser.parse_args()
    algorithm = AzimuthAlgorithm()
    result = algorithm.run(**args.__dict__)
    result.save()
