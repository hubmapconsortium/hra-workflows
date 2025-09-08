import argparse
import json
import typing as t

import anndata
import pandas as pd

from src.util.ids import create_cell_id


def get_unique_rows_with_counts(
    matrix: anndata.AnnData, clid_column: str
) -> pd.DataFrame:
    """Computes unique CLIDs and the total count for each.

    Args:
        matrix (anndata.AnnData): Data
        clid_column (str): Column with CLIDs

    Returns:
        pd.DataFrame: A frame with unique CLIDs and counts added
    """
    counts = matrix.obs.value_counts(clid_column).reset_index()
    counts.columns = [clid_column, "count"]
    obs_with_counts = matrix.obs.merge(counts, how="left")
    return obs_with_counts.drop_duplicates(clid_column)


def unique_rows_to_summary_rows(
    unique: pd.DataFrame,
    clid_column: str,
    label_column: str,
    gene_expr_column: str,
    nsforest_gene_expr_column: str,
    counts_column="count",
) -> t.List[dict]:
    """Converts a data frame with unique CLIDs rows into cell summary rows.

    Args:
        unique (pd.DataFrame): Data with unique CLIDs
        clid_column (str): Column with CLIDs
        label_column (str): Column with labels
        gene_expr_column (str): Column with gene expressions
        nsforest_gene_expr_column (str, optional): Column with NSForest gene expressions
        counts_column (str, optional): Column with the total counts. Defaults to "count".

    Returns:
        t.List[dict]: A cell summary for each row in the source data
    """
    columns = [clid_column, label_column, gene_expr_column, nsforest_gene_expr_column, counts_column]
    column_mapping = {
        clid_column: "cell_id",
        label_column: "cell_label",
        gene_expr_column: "gene_expr",
        nsforest_gene_expr_column: "nsforest_gene_expr",
        counts_column: "count",
    }
    
    df = unique[columns].rename(columns=column_mapping)

    df["@type"] = "CellSummaryRow"
    df["cell_id"] = df["cell_id"].map(create_cell_id)
    df["percentage"] = df["count"] / df["count"].sum()
    df["gene_expr"] = (
        df["gene_expr"]
        .astype(object)
        .apply(lambda x: [] if pd.isna(x) else json.loads(x))
    )
    df["nsforest_gene_expr"] = (
        df["nsforest_gene_expr"]
        .astype(object)
        .apply(lambda x: [] if pd.isna(x) else json.loads(x))
    ) 
    
    return df.to_dict("records")


def main(args: argparse.Namespace):
    """Extract and save a cell summary from annotated data.

    Args:
        args (argparse.Namespace):
            CLI arguments, must contain "matrix", "annotation_method",
            "cell_id_column", "cell_label_column", "gene_expr_column",
            "nsforest_gene_expr_column", "cell_source, "jsonld_context", "output", and
            "annotations_output"
    """
    unique_rows = get_unique_rows_with_counts(args.matrix, args.cell_id_column)
    summary_rows = unique_rows_to_summary_rows(
        unique_rows, 
        args.cell_id_column, 
        args.cell_label_column, 
        args.gene_expr_column,
        args.nsforest_gene_expr_column
    )
    summary = {
        "@type": "CellSummary",
        "annotation_method": args.annotation_method,
        "modality": "sc_transcriptomics",
        "cell_source": args.cell_source,
        "summary": summary_rows,
    }

    context: dict = json.load(args.jsonld_context)
    graph: list = context.setdefault("@graph", [])
    graph.append(summary)

    args.matrix.obs.to_csv(args.annotations_output, compression="gzip")

    
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
    parser.add_argument(
        "--nsforest-gene-expr-column",
        default="nsforest_gene_expr",
        help="NSForest gene expression column",
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
        type=str,
        default="annotations.csv.gz",
        help="Matrix obs",
    )

    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
