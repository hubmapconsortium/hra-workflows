import dataclasses
import enum
import io
import json
import pprint
import subprocess
import traceback
import typing as t
from pathlib import Path

import anndata


class Status(enum.Enum):
    SUCCESS = "success"
    FAILURE = "failure"


@dataclasses.dataclass
class AlgorithmReport:
    """Stores information about an algorithm run.

    Attributes:
        matrix (Path): Data h5ad output file path
        annotations (Path): Annotations csv file path
        report (Path): Report json file path
        status (Status): Algorithm run status
        data (anndata.AnnData, optional): Annotated data for a successful run
        failure_cause (t.Any): Data associated with a failed run
    """

    matrix: Path
    annotations: Path
    report: Path
    data: t.Optional[anndata.AnnData] = None
    status = Status.SUCCESS
    failure_cause: t.Any = None

    def is_success(self) -> bool:
        """Gets whether the run was successful.
        """
        return self.status == Status.SUCCESS

    def set_success(self, data: anndata.AnnData):
        """Set the run to success and add the result data to the report.

        Args:
            data (anndata.AnnData): Result data
        """
        self.status = Status.SUCCESS
        self.data = data
        self.failure_cause = None
        return self

    def set_failure(self, cause: t.Any):
        """Set the run to failure and add the cause to the report.

        Args:
            cause (t.Any): Any value but usually the error that caused the failure
        """
        self.status = Status.FAILURE
        self.data = anndata.AnnData()
        self.failure_cause = cause
        return self

    def save(self):
        """Saves the report and associated data to file.
        """
        self.save_matrix()
        self.save_report()

    def save_matrix(self):
        """Saves the matrix and annotations to file.
        """
        matrix = self.data
        if matrix is not None:
            matrix.obs.to_csv(self.annotations)
            matrix.write_h5ad(self.matrix)

    def save_report(self):
        """Saves the report status to file.
        """
        result = {"status": self.status.value}
        if not self.is_success():
            self.format_cause(result)

        with open(self.report, "w") as file:
            json.dump(result, file, indent=4)

    def format_cause(self, result: dict):
        """Format a failure cause and add it to a result dict.

        Args:
            result (dict): Dict to add formatted cause to
        """
        cause = self.failure_cause
        if isinstance(cause, Exception):
            result["cause"] = repr(cause)
            result["traceback"] = traceback.format_tb(cause.__traceback__)
            if isinstance(cause, subprocess.CalledProcessError):
                # stdout and stderr does not show in the repr of a CalledProcessError
                # Add them to the result here instead to give more context on errors
                result["stdout"] = cause.stdout
                result["stderr"] = cause.stderr
        else:
            stream = io.StringIO()
            pprint.pprint(cause, stream=stream)
            result["cause"] = stream.getvalue()
