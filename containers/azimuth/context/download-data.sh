#!/bin/bash
set -e

OUTPUT_DIR=${1:-"./azimuth"}
MAPPING_FILE=${2:-"/organ-mapping.json"}

mkdir -p "$OUTPUT_DIR"
Rscript /download_reference_data.R "$MAPPING_FILE" "$OUTPUT_DIR"
