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
        "gene_id": str,
        "gene_label": str,
        "count": int,
        "percentage": float,
    },
)

CellSummary = t.TypedDict(
    "CellSummary",
    {
        "@type": t.Literal["CellSummary"],
        "annotation_method": str,
        "cell_source": t.Optional[str],
        "summary": t.List[CellSummaryRow],
    },
)

_NON_WORDS_REGEX = re.compile("\W+")
_NON_ALHPA_NUM_HYPHEN_REGEX = re.compile("[^a-zA-Z0-9-]")


def _CellSummaryRow(
    cell_id: str, cell_label: str, gene_id: str, gene_label: str
) -> CellSummaryRow:
    return {
        "@type": "CellSummaryRow",
        "cell_id": cell_id,
        "cell_label": cell_label,
        "gene_id": gene_id,
        "gene_label": gene_label,
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


def _is_valid_column(reader: csv.DictReader, column: t.Optional[str]) -> bool:
    return column is None or column in reader.fieldnames


def normalize_id(id: str) -> str:
    id = id.lower().strip()
    id = _NON_WORDS_REGEX.sub("-", id)
    id = _NON_ALHPA_NUM_HYPHEN_REGEX.sub("", id)
    id = "ASCTB-TEMP:" + id
    return id


def compute_summary_rows(
    items: t.Iterator[t.Dict[str, str]],
    cell_id_column: t.Optional[str],
    cell_label_column: str,
    gene_id_column: t.Optional[str],
    gene_label_column: str,
) -> t.List[CellSummaryRow]:
    rowsById: t.Dict[t.Tuple[str, str], CellSummaryRow] = {}
    cell_id_transform = lambda val: val
    gene_id_transform = lambda val: val
    total = 0

    if cell_id_column is None:
        cell_id_column = cell_label_column
        cell_id_transform = normalize_id

    if gene_id_column is None:
        gene_id_column = gene_label_column
        gene_id_transform = normalize_id

    for item in items:
        id = (item[cell_id_column], item[gene_id_column])
        if id not in rowsById:
            rowsById[id] = _CellSummaryRow(
                cell_id_transform(item[cell_id_column]),
                item[cell_label_column],
                gene_id_transform(item[gene_id_column]),
                item[gene_label_column],
            )
        rowsById[id]["count"] += 1
        total += 1

    for row in rowsById.values():
        row["percentage"] = row["count"] / total

    return list(rowsById.values())


def main(args: argparse.Namespace):
    context: t.Dict = json.load(args.jsonld_context)
    reader = csv.DictReader(args.input)
    cell_id_column = args.cell_id_column
    cell_label_column = args.cell_label_column
    gene_id_column = args.gene_id_column
    gene_label_column = args.gene_label_column

    rows = []
    if (
        _is_valid_column(reader, cell_id_column)
        and _is_valid_column(reader, cell_label_column)
        and _is_valid_column(reader, gene_id_column)
        and _is_valid_column(reader, gene_label_column)
    ):
        rows = compute_summary_rows(
            reader, cell_id_column, cell_label_column, gene_id_column, gene_label_column
        )

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
        "--gene-label-column",
        required=True,
        help="Gene label column. Used for grouping if --gene-id-column is not provided.",
    )
    parser.add_argument(
        "--cell-id-column", help="Optional id column. Groups by label if not provided."
    )
    parser.add_argument(
        "--gene-id-column", help="Optional id column. Groups by label if not provided."
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
