import abc
import typing as t
from pathlib import Path

import anndata

from .organ import OrganLookup
from .report import AlgorithmReport

Organ = t.TypeVar("Organ")
Options = t.TypeVar("Options")
AnnDataOrPath = t.Union[str, Path, anndata.AnnData]
RunResult = t.Union[AnnDataOrPath, t.Tuple[AnnDataOrPath, str]]


class Algorithm(t.Generic[Organ, Options], abc.ABC):
    """An annotation algorithm.

    Attributes:
        organ_lookup (t.Callable[[Path], OrganLookup[Organ]]): Callable to create an organ lookup
        prediction_column (t.Optional[str]): Column in annotated data with the predictions
    """

    def __init__(
        self,
        organ_lookup: t.Callable[[Path], OrganLookup[Organ]],
        prediction_column: t.Optional[str] = None,
    ):
        self.organ_lookup = organ_lookup
        self.prediction_column = prediction_column

    def run(
        self,
        matrix: Path,
        organ: str,
        organ_mapping: Path,
        output_matrix: Path,
        output_annotations: Path,
        output_report: Path,
        **options,
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
            lookup = self.organ_lookup(organ_mapping)
            result = self.do_run(matrix, lookup.get(organ), t.cast(Options, options))
            result = self.__post_process_result(result)
            report.set_success(result)
        except Exception as error:
            report.set_failure(error)

        return report

    @abc.abstractmethod
    def do_run(self, matrix: Path, organ: Organ, options: Options) -> RunResult:
        """Perform a annotation run. Must be overridden in subclasses.

        Args:
            matrix (Path): Path to the h5ad data file
            organ (Organ): Organ associated with the data
            options (Options): Additional algorithm specific options

        Returns:
            RunResult:
                Annotated data either in-memory or a path to a h5ad.
                Can also return a tuple where the first element is
                the annotated data and the second element is the name
                of the column that stores the predictions.
        """
        ...

    def __post_process_result(self, result: RunResult) -> anndata.AnnData:
        """Normalize the result of a run.

        Args:
            result (RunResult): Non-normalized result value

        Returns:
            anndata.AnnData: Loaded h5ad data
        """
        prediction_column = self.prediction_column
        if isinstance(result, tuple):
            result, prediction_column = result
        if isinstance(result, (str, Path)):
            result = anndata.read_h5ad(result)
        if prediction_column is not None and prediction_column in result.obs.columns:
            result.obs["hra_prediction"] = result.obs[prediction_column]

        return result
