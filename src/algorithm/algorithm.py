import abc
import json
import typing as t
from pathlib import Path

import anndata

from .report import AlgorithmReport

OrganMetadata = t.TypeVar("OrganMetadata", bound=dict)
Options = t.TypeVar("Options", bound=dict)
AnnDataOrPath = t.Union[str, Path, anndata.AnnData]


class RunResult(t.TypedDict):
    data: AnnDataOrPath
    organ_level: str
    prediction_column: t.Optional[str]


class Algorithm(t.Generic[OrganMetadata, Options], abc.ABC):
    """An annotation algorithm.

    Attributes:
        prediction_column (t.Optional[str]): Column in annotated data with the predictions
    """

    def __init__(self, prediction_column: t.Optional[str] = None):
        self.prediction_column = prediction_column

    def run(
        self,
        matrix: Path,
        organ: str,
        organ_metadata: Path,
        output_matrix: Path,
        output_annotations: Path,
        output_report: Path,
        **options: Options,
    ) -> AlgorithmReport:
        """Runs the algorithm to annotate data.

        Args:
            matrix (Path): Path to h5ad data file
            organ (str): Raw organ identifier
            organ_mapping (Path): Path to json file containing organ mapping information
            output_matrix (Path): Path where the annotated h5ad file will be written
            output_annotations (Path): Path where the annotation csv will be written
            output_report (Path): Path where the algorithm report json will be written

        Returns:
            AlgorithmReport: Report containing the status of the run
        """
        report = AlgorithmReport(output_matrix, output_annotations, output_report)
        try:
            resolved_organ, metadata = self.__load_metadata(organ, organ_metadata)
            result = self.do_run(matrix, organ, metadata, t.cast(Options, options))
            data = self.__post_process_result(result, resolved_organ, metadata)
            report.set_success(data)
        except Exception as error:
            report.set_failure(error)

        return report

    @abc.abstractmethod
    def do_run(
        self, matrix: Path, organ: str, metadata: OrganMetadata, options: Options
    ) -> RunResult:
        """Perform a annotation run. Must be overridden in subclasses.

        Args:
            matrix (Path): Path to the h5ad data file
            organ (Organ): Organ associated with the data
            options (Options): Additional algorithm specific options

        Returns:
            RunResult:
                Annotated data either in-memory or a path to a h5ad,
                the organ level to use in crosswalking, and
                optionally the column storing predictions
        """
        ...

    def __load_metadata(
        self, organ: str, organ_metadata: Path
    ) -> t.Tuple[str, OrganMetadata]:
        """Loads metadata for an organ from file.

        Args:
            organ (str): Organ id
            organ_metadata (Path): Path to metadata file

        Returns:
            t.Tuple[str, OrganMetadata]: The resolved organ id and associated metadata

        Raises:
            ValueError: If the organ is not supported by the algorithm
        """
        with open(organ_metadata) as file:
            metadata = json.load(file)
        return self.__resolve_metadata(organ, organ, metadata, set())

    def __resolve_metadata(
        self,
        original_organ: str,
        organ: str,
        metadata: t.Dict[str, t.Union[str, OrganMetadata]],
        seen: t.Set[str],
    ) -> t.Tuple[str, OrganMetadata]:
        """Resolves organ and metadata following links in the metadata table.

        Args:
            original_organ (str): Original organ id
            organ (str): Organ id
            metadata (t.Dict[str, t.Union[str, OrganMetadata]]): Loaded metadata table
            seen (t.Set[str]): Set of seen organ ids used to detect cycles in the metadata table

        Returns:
            t.Tuple[str, OrganMetadata]: The resolved organ id and associated metadata

        Raises:
            ValueError: If resolving metadata for the organ fails
        """
        if organ not in metadata:
            raise ValueError(f"Organ {original_organ} is not supported")
        if organ in seen:
            raise ValueError(
                f"Circular metadata links detected for organ {original_organ}"
            )

        seen.add(organ)
        value = metadata[organ]
        if isinstance(value, str):
            return self.__resolve_metadata(original_organ, value, metadata, seen)
        return organ, value

    def __post_process_result(
        self, result: RunResult, organ: str, metadata: OrganMetadata
    ) -> anndata.AnnData:
        """Normalize the result of a run.

        Args:
            result (RunResult): Run result dictionary
            organ (str): Organ id
            metadata (str): Organ metadata

        Returns:
            anndata.AnnData: Loaded h5ad data
        """

        data = result["data"]
        if isinstance(data, (str, Path)):
            data = anndata.read_h5ad(data)

        prediction_column = result.get("prediction_column", self.prediction_column)
        if prediction_column is None:
            raise ValueError("Missing prediction column")
        elif prediction_column not in data.obs.columns:
            raise ValueError(
                "Prediction column does not exist in the result", prediction_column
            )
        else:
            data.obs["hra_prediction"] = data.obs[prediction_column]

        data.uns["hra_organ_metadata"] = metadata
        data.uns["hra_crosswalking"] = {
            "organ_id": organ,
            "organ_level": result["organ_level"],
        }

        return data
