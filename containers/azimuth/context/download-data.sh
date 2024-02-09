#!/bin/bash
set -e

OUTPUT_DIR=${1:-"./azimuth"}
METADATA_FILE=${2:-"/organ-metadata.json"}
export R_LIBS="$OUTPUT_DIR"

mkdir -p "$OUTPUT_DIR"
Rscript /download_reference_data.R "$METADATA_FILE" "$OUTPUT_DIR"
