import argparse
import json

import anndata
import pandas as pd


def get_unique_rows_with_counts(
    matrix: anndata.AnnData, clid_column: str
) -> pd.DataFrame:
    counts = matrix.obs.value_counts(clid_column).reset_index()
    counts.columns = [clid_column, "count"]
    obs_with_counts = matrix.obs.merge(counts, how="left")
    return obs_with_counts.drop_duplicates(clid_column)


def unique_rows_to_summary_rows(
    unique: pd.DataFrame,
    clid_column: str,
    label_column: str,
    gene_expr_column: str,
    counts_column="count",
):
    columns = [clid_column, label_column, gene_expr_column, counts_column]
    df = unique[columns].rename(
        columns={
            clid_column: "cell_id",
            label_column: "cell_label",
            gene_expr_column: "gene_expr",
            counts_column: "count",
        }
    )

    df["@type"] = "CellSummaryRow"
    df["percentage"] = df["count"] / df["count"].sum()
    df["gene_expr"] = df["gene_expr"].astype(object).apply(lambda x: [] if pd.isna(x) else json.loads(x))
    return df.to_dict("records")


def main(args: argparse.Namespace):
    unique_rows = get_unique_rows_with_counts(args.matrix, args.cell_id_column)
    summary_rows = unique_rows_to_summary_rows(
        unique_rows, args.cell_id_column, args.cell_label_column, args.gene_expr_column
    )
    summary = {
        "@type": "CellSummary",
        "annotation_method": args.annotation_method,
        "modality": "bulk",
        "cell_source": args.cell_source,
        "summary": summary_rows,
    }

    context: dict = json.load(args.jsonld_context)
    graph: list = context.setdefault("@graph", [])
    graph.append(summary)

    args.matrix.obs.to_csv(args.annotations_output)
    json.dump(context, args.output, indent=2)


def _get_arg_parser():
    parser = argparse.ArgumentParser(
        description="Extract cell summary from annotations"
    )
    parser.add_argument("matrix", type=anndata.read_h5ad, help="Annotated h5ad matrix")
    parser.add_argument(
        "--annotation-method", required=True, help="Method used to extract annotations"
    )
    parser.add_argument(
        "--cell-id-column",
        default="clid",
        help="Cell id column",
    )
    parser.add_argument(
        "--cell-label-column",
        default="hra_prediction",
        help="Cell label column",
    )
    parser.add_argument(
        "--gene-expr-column",
        default="gene_expr",
        help="Gene expression column",
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
    parser.add_argument(
        "--annotations-output",
        type=argparse.FileType("w"),
        default="annotations.csv",
        help="Matrix obs",
    )

    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
