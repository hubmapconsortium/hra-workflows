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
        ...

    def __post_process_result(self, result: RunResult) -> anndata.AnnData:
        prediction_column = self.prediction_column
        if isinstance(result, tuple):
            result, prediction_column = result
        if isinstance(result, (str, Path)):
            result = anndata.read_h5ad(result)
        if prediction_column is not None and prediction_column in result.obs.columns:
            result.obs["hra_prediction"] = result.obs[prediction_column]

        return result
