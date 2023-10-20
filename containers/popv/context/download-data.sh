#!/bin/bash
set -e

OUTPUT_DIR=${3:-"./popv"}

mkdir -p "$OUTPUT_DIR"
/download-models.sh "$1" "$OUTPUT_DIR/models"
/download-reference-data.sh "$2" "$OUTPUT_DIR/reference-data"
