import traceback
import dataclasses
import enum
import io
import json
import pprint
import typing as t
from pathlib import Path

import anndata
import pandas


class Status(enum.Enum):
    SUCCESS = "success"
    FAILURE = "failure"


@dataclasses.dataclass
class AlgorithmReport:
    matrix: Path
    annotations: Path
    report: Path
    data: t.Union[Path, anndata.AnnData]
    status = Status.SUCCESS
    failure_cause: t.Any = None

    def is_success(self) -> bool:
        return self.status == Status.SUCCESS

    def set_success(self, data: t.Union[Path, anndata.AnnData]):
        self.status = Status.SUCCESS
        self.data = data
        self.failure_cause = None
        return self

    def set_failure(self, cause: t.Any):
        self.status = Status.FAILURE
        self.failure_cause = cause
        return self

    def save(self):
        self.save_matrix()
        self.save_report()

    def save_matrix(self):
        if isinstance(self.data, Path):
            self.data = anndata.read_h5ad(self.data)
        self.data.write_h5ad(self.matrix)

        obs = self.data.obs
        if not self.is_success():
            # Create an empty observations frame with the same columns as the original
            obs = pandas.DataFrame(columns=obs.columns)
        obs.to_csv(self.annotations)

    def save_report(self):
        result = {"status": self.status.value}
        if not self.is_success():
            self.format_cause(result)

        with open(self.report, 'w') as file:
            json.dump(result, file, indent=4)

    def format_cause(self, result: dict):
        cause = self.failure_cause
        if isinstance(cause, Exception):
            result["cause"] = repr(cause)
            result["traceback"] = traceback.format_tb(cause.__traceback__)
        else:
            stream = io.StringIO()
            pprint.pprint(cause, stream=stream)
            result["cause"] = stream.getvalue()
