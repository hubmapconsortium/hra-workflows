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
    gene_expr_data: t.Optional[dict] = None,
    nsforest_gene_expr_data: t.Optional[dict] = None,
    counts_column="count",
) -> t.List[dict]:
    """Converts a data frame with unique CLIDs rows into cell summary rows.

    Args:
        unique (pd.DataFrame): Data with unique CLIDs
        clid_column (str): Column with CLIDs
        label_column (str): Column with labels
        gene_expr_data (dict, optional): Gene expression data from JSON file
        nsforest_gene_expr_data (dict, optional): NSForest gene expression data from JSON file
        counts_column (str, optional): Column with the total counts. Defaults to "count".

    Returns:
        t.List[dict]: A cell summary for each row in the source data
    """
    # Get basic info from matrix
    columns = [clid_column, label_column, counts_column]
    column_mapping = {
        clid_column: "cell_id",
        label_column: "cell_label",
        counts_column: "count",
    }

    df = unique[columns].rename(columns=column_mapping)

    # Add standard fields
    df["@type"] = "CellSummaryRow"
    df["cell_id"] = df["cell_id"].map(create_cell_id)
    df["percentage"] = df["count"] / df["count"].sum()

    # Add gene expression data from JSON if provided
    if gene_expr_data is not None:
        df["gene_expr"] = unique[clid_column].astype(object).apply(lambda x: gene_expr_data.get(x, []))

    # Add NSForest gene expression data from JSON if provided
    if nsforest_gene_expr_data is not None:
        df["nsforest_gene_expr"] = unique[clid_column].astype(object).apply(lambda x: nsforest_gene_expr_data.get(x, []))

    return df.to_dict("records")


def main(args: argparse.Namespace):
    """Extract and save a cell summary from annotated data.

    Args:
        args (argparse.Namespace):
            CLI arguments, must contain "matrix", "annotation_method",
            "cell_id_column", "cell_label_column", "gene_expr_json",
            "nsforest_gene_expr_json", "cell_source", "jsonld_context", "output"
    """
    # Load gene expression data from JSON files
    gene_expr_data = None
    nsforest_gene_expr_data = None
    
    if args.gene_expr_json:
        with open(args.gene_expr_json) as f:
            gene_expr_data = json.load(f)
    if args.nsforest_gene_expr_json:
        with open(args.nsforest_gene_expr_json) as f:
            nsforest_gene_expr_data = json.load(f)
            
    # Get unique rows with counts from matrix
    unique_rows = get_unique_rows_with_counts(args.matrix, args.cell_id_column)
    
    # Create summary rows using JSON data
    summary_rows = unique_rows_to_summary_rows(
        unique_rows, 
        args.cell_id_column, 
        args.cell_label_column,
        gene_expr_data=gene_expr_data,
        nsforest_gene_expr_data=nsforest_gene_expr_data
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

    json.dump(context, args.output, indent=2)

def _get_arg_parser():
    parser = argparse.ArgumentParser(
        description="Extract cell summary from annotations"
    )
    parser.add_argument("matrix", type=anndata.read_h5ad, help="Matrix h5ad for cell labels and counts")
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
        "--gene-expr-json",
        type=str,
        help="Gene expression data JSON file",
    )
    parser.add_argument(
        "--nsforest-gene-expr-json",
        type=str,
        help="NSForest gene expression data JSON file",
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
