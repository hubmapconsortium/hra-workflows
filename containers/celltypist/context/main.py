from argparse import ArgumentParser
from logging import captureWarnings, getLogger
from logging.config import fileConfig
from pathlib import Path
from typing import Optional

import celltypist
from scanpy import AnnData
from scanpy import pp as preprocess
from scanpy import read_h5ad

_logger = getLogger(__name__)

def _get_arg_parser():
    parser = ArgumentParser(description='Compute annotations using celltypist')
    parser.add_argument('filename', type=Path, help='h5ad data file')
    parser.add_argument('--model', required=True, help='Model name. See https://www.celltypist.org/models')
    parser.add_argument('--existing-annotations-column', help='Column with existing annotations')
    parser.add_argument('--out', type=Path, help='Output file', default='annotations.csv')
    parser.add_argument('--logger-config', type=Path, help='Python logging configuration file')

    return parser

def normalize_data(data: AnnData):
    primary_column = 'feature_name'
    alternative_primary_column = 'gene_symbol'
    if primary_column not in data.var.columns:
        _logger.warning(f'Missing {primary_column} data column. Attempting to use {alternative_primary_column} instead')
        data.var = data.var.rename({ alternative_primary_column: primary_column })

    # TODO I've excluded some of the original code from an if statement that is always false
    # https://github.com/cns-iu/ct-ann-predictive-analytics/blob/main/celltypist/celltypist_prediction_pipeline.py#L116

    preprocess.normalize_total(data, target_sum=1e4)
    preprocess.log1p(data)
    data.var_names_make_unique()

def main(filename: Path, model_name: str, existing_annotations_column: Optional[str], output_filename: Path):
    data = read_h5ad(filename)
    model = celltypist.models.Model.load(model_name)

    normalize_data(data)
    predictions = celltypist.annotate(data, model, majority_voting=True)
    annotations = predictions.to_adata()

    if existing_annotations_column:
        matches = lambda row: row[existing_annotations_column] == row['majority_voting']
        annotations.obs['exp_vs_pred'] = annotations.obs.apply(matches, axis=1)

    annotations.obs.to_csv(output_filename, index=True)

if __name__ == '__main__':
    parser = _get_arg_parser()
    args = parser.parse_args()
    if args.logger_config is not None:
        captureWarnings(True)
        fileConfig(args.logger_config, disable_existing_loggers=False)
    main(args.filename, args.model, args.existing_annotations_column, args.out)
