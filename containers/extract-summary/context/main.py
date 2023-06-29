import argparse
import csv
import json
import re
import typing as t

CellSummaryRow = t.TypedDict(
    "CellSummaryRow",
    {
        "@type": t.Literal["CellSummaryRow"],
        "cell_id": str,
        "cell_label": str,
        "count": int,
        "percentage": float,
    },
)

CellSummary = t.TypedDict(
    "CellSummary",
    {
        "@type": t.Literal["CellSummary"],
        "annotation_method": str,
        "cell_source": t.NotRequired[str],
        "summary": t.List[CellSummaryRow],
    },
)

_NON_WORDS_REGEX = re.compile("\W+")
_NON_ALHPA_NUM_HYPHEN_REGEX = re.compile("[^a-zA-Z0-9-]")


def _CellSummaryRow(id: str, label: str) -> CellSummaryRow:
    return {
        "@type": "CellSummaryRow",
        "cell_id": id,
        "cell_label": label,
        "count": 0,
        "percentage": 0,
    }


def _CellSummary(
    method: str, source: t.Optional[str], summaries: t.List[CellSummaryRow]
) -> CellSummary:
    result: CellSummary = {
        "@type": "CellSummary",
        "annotation_method": method,
        "cell_source": source,
        "summary": summaries,
    }

    if source is None:
        del result["cell_source"]
    return result


def _check_column(reader: csv.DictReader, column: t.Optional[str]):
    if column is not None and column not in reader.fieldnames:
        raise ValueError(f"{column} is not a column in the csv file")


def normalize_id(id: str) -> str:
    id = id.lower().strip()
    id = _NON_WORDS_REGEX.sub("-", id)
    id = _NON_ALHPA_NUM_HYPHEN_REGEX.sub("", id)
    id = "ASCTB-TEMP:" + id
    return id


def compute_summary_rows(
    items: t.Iterator[t.Dict[str, str]], id_column: t.Optional[str], label_column: str
) -> t.List[CellSummaryRow]:
    rowsById: t.Dict[str, CellSummaryRow] = {}
    id_transform = lambda val: val
    total = 0

    if id_column is None:
        id_column = label_column
        id_transform = normalize_id

    for item in items:
        id = item[id_column]
        if id not in rowsById:
            rowsById[id] = _CellSummaryRow(id_transform(id), item[label_column])
        rowsById[id]["count"] += 1
        total += 1

    for row in rowsById.values():
        row["percentage"] = row["count"] / total

    return list(rowsById.values())


def main(args: argparse.Namespace):
    context: t.Dict = json.load(args.jsonld_context)
    reader = csv.DictReader(args.input)
    id_column = args.cell_id_column
    label_column = args.cell_label_column

    _check_column(reader, id_column)
    _check_column(reader, label_column)

    rows = compute_summary_rows(reader, id_column, label_column)
    summary = _CellSummary(args.annotation_method, args.cell_source, rows)
    context.setdefault("@graph", []).append(summary)
    json.dump(context, args.output, indent=2)


def _get_arg_parser():
    parser = argparse.ArgumentParser(
        description="Extract cell summary from annotations"
    )
    parser.add_argument(
        "input", type=argparse.FileType("r"), help="Annotations csv file"
    )
    parser.add_argument(
        "--annotation-method", required=True, help="Method used to extract annotations"
    )
    parser.add_argument(
        "--cell-label-column",
        required=True,
        help="Cell label column. Used for grouping if --cell-id-column is not provided.",
    )
    parser.add_argument(
        "--cell-id-column", help="Optional id column. Groups by label if not provided."
    )
    parser.add_argument("--cell-source", help="Cell source. Must be an IRI.")
    parser.add_argument(
        "--jsonld-context",
        type=argparse.FileType("r"),
        default="/context.jsonld",
        help="Base jsonld context to add summary to.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=argparse.FileType("w"),
        default="summary.jsonld",
        help="Output file",
    )

    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
