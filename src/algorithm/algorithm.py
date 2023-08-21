import abc
import typing as t
from pathlib import Path

import anndata

from .organ import OrganLookup
from .report import AlgorithmReport

Organ = t.TypeVar("Organ")
Options = t.TypeVar("Options")


class Algorithm(t.Generic[Organ, Options], abc.ABC):
    def __init__(self, organ_lookup: t.Callable[[Path], OrganLookup[Organ]]):
        self.organ_lookup = organ_lookup

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
        report = AlgorithmReport(
            output_matrix, output_annotations, output_report, matrix
        )
        try:
            lookup = self.organ_lookup(organ_mapping)
            result = self.do_run(matrix, lookup.get(organ), t.cast(Options, options))
            report.set_success(result)
        except Exception as error:
            report.set_failure(error)

        return report

    @abc.abstractmethod
    def do_run(
        self, matrix: Path, organ: Organ, options: Options
    ) -> t.Union[Path, anndata.AnnData]:
        ...
