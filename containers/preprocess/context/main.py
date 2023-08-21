import argparse
from pathlib import Path

import anndata


def add_genes(
    matrix: anndata.AnnData, gene_column: str, gene_output_column: str
) -> anndata.AnnData:
    assert matrix.X is not None  # Makes type checker happy on the next line :)
    indices = matrix.X.argmax(axis=1)
    genes = matrix.var[gene_column][indices.flat]
    matrix.obs[gene_output_column] = genes.values
    return matrix


def main(args: argparse.Namespace):
    matrix = add_genes(args.matrix, args.gene_column, args.gene_output_column)
    matrix.write_h5ad(args.output_matrix)


def _get_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Preprocess h5ad data")
    parser.add_argument("matrix", type=anndata.read_h5ad, help="h5ad data file")
    parser.add_argument("--gene-column", help="Column with gene names")
    parser.add_argument(
        "--gene-output-column", default="gene", help="Output column for gene"
    )
    parser.add_argument(
        "--output-matrix",
        type=Path,
        default="preprocessed_matrix.h5ad",
        help="Preprocessed matrix output path",
    )

    return parser


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
