#!/usr/bin/env python3
"""
qc_metrics.py

Compute standard single-cell RNA-seq quality control (QC) metrics from an .h5ad file
using AnnData and Scanpy.

Features:
- Computes mitochondrial and ribosomal gene percentages
- Flags low-quality cells based on user-defined thresholds
- Exports per-cell and summary QC tables
- Optionally filters out low-quality cells and writes a filtered .h5ad
- Writes a machine-readable JSON summary for pipelines
- Usable both as a CLI and an importable Python module
- Uses Python logging (numeric levels, default INFO=20)
"""

from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from src.util.ensemble import add_ensemble_data


# -----------------------------------------------------------------------------
# Gene list loading
# -----------------------------------------------------------------------------
def _load_gene_list_from_csv(csv_path: str) -> Tuple[Set[str], Set[str]]:
    """
    Load gene names and Ensembl IDs from a CSV file.
    
    Parameters
    ----------
    csv_path : str
        Path to the CSV file containing gene information.
    
    Returns
    -------
    Tuple[Set[str], Set[str]]
        A tuple of (gene_names, ensembl_ids_stripped) where ensembl_ids_stripped
        has the version removed (everything after the dot).
    """
    try:
        df = pd.read_csv(csv_path)
        gene_names = set(df["gene_name"].dropna().unique())
        # Strip version from Ensembl IDs (remove everything after the dot)
        ensembl_ids = set(
            eid.split(".")[0] for eid in df["ensembl_gene_id"].dropna().unique()
        )
        logging.info(f"Loaded {len(gene_names)} genes and {len(ensembl_ids)} Ensembl IDs from {csv_path}")
        return gene_names, ensembl_ids
    except FileNotFoundError:
        logging.error(f"Gene list file not found: {csv_path}")
        raise
    except Exception as e:
        logging.error(f"Error loading gene list from {csv_path}: {e}")
        raise


# -----------------------------------------------------------------------------
# Logging setup
# -----------------------------------------------------------------------------
def _setup_logging(level: int = 20) -> None:
    """
    Configure logging format and level.

    Parameters
    ----------
    level : int
        Logging level (10=DEBUG, 20=INFO, 30=WARNING, 40=ERROR, 50=CRITICAL).
    """
    logging.basicConfig(
        format="%(levelname)s: %(message)s",
        level=level,
        datefmt="%Y-%m-%d %H:%M:%S",
    )


# -----------------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------------
def _get_arg_parser() -> argparse.ArgumentParser:
    """Return a configured argument parser for QC metrics script."""
    # Get default CSV paths from the same directory as this script
    script_dir = Path(__file__).parent
    default_mt_csv = script_dir / "human_mitochondrial_genes.csv"
    default_ribo_csv = script_dir / "human_ribosomal_genes_structural.csv"
    
    parser = argparse.ArgumentParser(
        description="Compute QC metrics from an h5ad file and optionally filter low-quality cells."
    )
    parser.add_argument("input", help="Path to the input .h5ad file.")
    parser.add_argument(
        "-o", "--output", help="Output directory for QC results.", default="qc_results"
    )
    parser.add_argument(
        "--mt-csv",
        help="Path to CSV file with mitochondrial gene information (gene_name and ensembl_gene_id columns required).",
        default=str(default_mt_csv),
    )
    parser.add_argument(
        "--ribo-csv",
        help="Path to CSV file with ribosomal gene information (gene_name and ensembl_gene_id columns required).",
        default=str(default_ribo_csv),
    )
    parser.add_argument(
        "--filter",
        action="store_true",
        help="Filter low-quality cells and save a filtered .h5ad file.",
    )
    parser.add_argument(
        "--log-level",
        type=int,
        default=20,
        help="Set the logging level numerically (10=DEBUG, 20=INFO, 30=WARNING, 40=ERROR, 50=CRITICAL). Default: 20 (INFO).",
    )
    return parser


# -----------------------------------------------------------------------------
# QC metric computation
# -----------------------------------------------------------------------------
def compute_qc_metrics(
    adata: ad.AnnData, 
    mt_genes: Optional[Tuple[Set[str], Set[str]]] = None,
    ribo_genes: Optional[Tuple[Set[str], Set[str]]] = None,
) -> ad.AnnData:
    """
    Compute per-cell QC metrics for mitochondrial and ribosomal genes.
    
    Parameters
    ----------
    adata : ad.AnnData
        Annotated data matrix.
    mt_genes : Optional[Tuple[Set[str], Set[str]]]
        Tuple of (gene_names, ensembl_ids_stripped) for mitochondrial genes.
        If None, no mitochondrial genes will be marked.
    ribo_genes : Optional[Tuple[Set[str], Set[str]]]
        Tuple of (gene_names, ensembl_ids_stripped) for ribosomal genes.
        If None, no ribosomal genes will be marked.
    
    Returns
    -------
    ad.AnnData
        Updated AnnData object with QC metrics.
    """
    def _is_in_gene_set(var_name: str, gene_set: Optional[Tuple[Set[str], Set[str]]]) -> bool:
        """Check if a variable name matches a gene name or Ensembl ID (version-stripped)."""
        if gene_set is None:
            return False
        gene_names, ensembl_ids = gene_set
        # Check exact match on gene name (case-sensitive)
        if var_name in gene_names:
            return True
        # Check if it's an Ensembl ID (strip version for comparison)
        var_name_stripped = var_name.split(".")[0]
        if var_name_stripped in ensembl_ids:
            return True
        return False

    # Mark genes based on gene lists
    adata.var["mt"] = adata.var_names.map(lambda x: _is_in_gene_set(x, mt_genes))
    adata.var["ribo"] = adata.var_names.map(lambda x: _is_in_gene_set(x, ribo_genes))

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"], percent_top=(20, 50), log1p=False, inplace=True
    )
    return adata


def flag_low_quality_cells(
    adata: ad.AnnData,
    min_genes: int,
    max_genes: int,
    max_mt: float,
    min_ribo: Optional[float] = None,
) -> pd.DataFrame:
    """Flag low-quality cells based on QC thresholds."""
    low_quality = (
        (adata.obs["n_genes_by_counts"] < min_genes)
        | (adata.obs["n_genes_by_counts"] > max_genes)
        | (adata.obs["pct_counts_mt"] > max_mt)
    )

    if min_ribo is not None:
        low_quality |= adata.obs["pct_counts_ribo"] < min_ribo

    adata.obs["low_quality"] = low_quality

    qc_df = adata.obs[
        [
            "n_genes_by_counts",
            "total_counts",
            "pct_counts_in_top_20_genes",
            "pct_counts_in_top_50_genes",
            "pct_counts_mt",
            "pct_counts_ribo",
            "low_quality",
        ]
    ].copy()

    return qc_df


# -----------------------------------------------------------------------------
# Reporting and persistence
# -----------------------------------------------------------------------------
def summarize_qc(qc_df: pd.DataFrame) -> pd.DataFrame:
    """Generate summary statistics for QC metrics."""
    summary = qc_df.describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95]).T
    summary.loc["mean_n_genes_by_counts", "mean"] = qc_df["n_genes_by_counts"].mean()
    summary.loc["mean_total_counts", "mean"] = qc_df["total_counts"].mean()
    summary.loc["mean_pct_counts_in_top_20_genes", "mean"] = qc_df["pct_counts_in_top_20_genes"].mean()
    summary.loc["mean_pct_counts_in_top_50_genes", "mean"] = qc_df["pct_counts_in_top_50_genes"].mean()
    summary.loc["mean_pct_counts_mt", "mean"] = qc_df["pct_counts_mt"].mean()
    summary.loc["mean_pct_counts_ribo", "mean"] = qc_df["pct_counts_ribo"].mean()
    return summary


def write_qc_outputs(
    qc_df: pd.DataFrame, summary: pd.DataFrame, output_dir: str
) -> Tuple[str, str]:
    """Write per-cell and summary QC results to CSV files."""
    os.makedirs(output_dir, exist_ok=True)
    per_cell_path = os.path.join(output_dir, "qc_per_cell.csv")
    summary_path = os.path.join(output_dir, "qc_summary.csv")

    qc_df.to_csv(per_cell_path)
    summary.to_csv(summary_path)

    logging.info(f"Wrote QC per-cell metrics: {per_cell_path}")
    logging.info(f"Wrote QC summary statistics: {summary_path}")

    return per_cell_path, summary_path


def filter_and_save(adata: ad.AnnData, output_dir: str) -> str:
    """Filter out low-quality cells and write a filtered .h5ad file."""
    filtered_path = os.path.join(output_dir, "filtered_data.h5ad")
    adata[~adata.obs["low_quality"]].write(filtered_path)
    logging.info(f"Filtered AnnData file written: {filtered_path}")
    return filtered_path


def write_json_summary(
    output_dir: str, metadata: Dict[str, Any], file_paths: Dict[str, str]
) -> str:
    """Write a JSON summary report with key QC statistics and file paths."""

    def _convert(o):
        if isinstance(o, (np.generic,)):
            return o.item()
        elif isinstance(o, dict):
            return {k: _convert(v) for k, v in o.items()}
        elif isinstance(o, list):
            return [_convert(v) for v in o]
        return o

    json_path = os.path.join(output_dir, "qc_summary.json")
    summary_data = _convert({**metadata, "files": file_paths})

    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(summary_data, f, indent=4)

    logging.info(f"JSON summary written: {json_path}")
    return json_path


# -----------------------------------------------------------------------------
# Main workflow (importable)
# -----------------------------------------------------------------------------
def run_qc(
    input_path: str,
    output_dir: str = "qc_results",
    mt_csv: Optional[str] = None,
    ribo_csv: Optional[str] = None,
    filter_cells: bool = False,
) -> Dict[str, str]:
    """Run the QC pipeline programmatically."""
    logging.info(f"Loading data from {input_path}")
    adata = sc.read_h5ad(input_path)

    # Load gene lists from CSV files if provided
    mt_genes = None
    ribo_genes = None
    
    if mt_csv:
        logging.info(f"Loading mitochondrial genes from {mt_csv}")
        mt_genes = _load_gene_list_from_csv(mt_csv)
    
    if ribo_csv:
        logging.info(f"Loading ribosomal genes from {ribo_csv}")
        ribo_genes = _load_gene_list_from_csv(ribo_csv)

    logging.info("Computing QC metrics ...")
    adata = compute_qc_metrics(adata, mt_genes=mt_genes, ribo_genes=ribo_genes)

    logging.info("Flagging low-quality cells ...")
    qc_df = flag_low_quality_cells(
        adata, 
        min_genes=200,
        max_genes=7500,
        max_mt=5.0,
        min_ribo=None,
    )
    summary = summarize_qc(qc_df)

    per_cell_path, summary_path = write_qc_outputs(qc_df, summary, output_dir)

    total_cells = adata.n_obs
    low_quality_cells = int(adata.obs["low_quality"].sum())
    pct_removed = (low_quality_cells / total_cells) * 100
    mean_mt = float(summary.loc['pct_counts_mt', 'mean'])
    mean_ribo = float(summary.loc['pct_counts_ribo', 'mean'])

    logging.info(f"Total cells: {total_cells}")
    logging.info(f"Flagged low-quality cells: {low_quality_cells} ({pct_removed:.2f}%)")
    logging.info(f"Mean mitochondrial %: {mean_mt:.2f}")
    logging.info(f"Mean ribosomal %: {mean_ribo:.2f}")

    results = {
        "qc_per_cell_csv": per_cell_path,
        "qc_summary_csv": summary_path,
    }

    if filter_cells:
        logging.info("Filtering low-quality cells ...")
        filtered_path = filter_and_save(adata, output_dir)
        results["filtered_h5ad"] = filtered_path
    else:
        logging.info("Filtering skipped (use filter_cells=True or --filter).")

    metadata = {
        "input_file": os.path.abspath(input_path),
        "total_cells": int(total_cells),
        "low_quality_cells": int(low_quality_cells),
        "percent_low_quality": round(pct_removed, 2),
        "mean_n_genes_by_counts": round(summary.loc['n_genes_by_counts', 'mean'], 3),
        "mean_total_counts": round(summary.loc['total_counts', 'mean'], 3),
        "mean_pct_counts_in_top_20_genes": round(summary.loc['pct_counts_in_top_20_genes', 'mean'], 3),
        "mean_pct_counts_in_top_50_genes": round(summary.loc['pct_counts_in_top_50_genes', 'mean'], 3),
        "mean_pct_counts_mt": round(mean_mt, 3),
        "mean_pct_counts_ribo": round(mean_ribo, 3),
        "thresholds": {
            "min_genes": 200,
            "max_genes": 7500,
            "max_mt": 5.0,
            "min_ribo": None,
            "mt_csv": mt_csv,
            "ribo_csv": ribo_csv,
        },
    }

    json_path = write_json_summary(output_dir, metadata, results)
    results["qc_summary_json"] = json_path

    logging.info("QC analysis complete.")
    return results


# -----------------------------------------------------------------------------
# CLI entry point
# -----------------------------------------------------------------------------
def main() -> None:
    """Entry point for command-line usage."""
    parser = _get_arg_parser()
    args = parser.parse_args()
    _setup_logging(args.log_level)

    run_qc(
        input_path=args.input,
        output_dir=args.output,
        mt_csv=args.mt_csv,
        ribo_csv=args.ribo_csv,
        filter_cells=args.filter,
    )


if __name__ == "__main__":
    main()
