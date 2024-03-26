#!/bin/bash
set -e

REFERENCE_DATA_ID=${1:?"A zenodo reference data id must be provided to download!"}
REFERENCE_DATA_DIR=${2:-"./popv/reference-data"}

mkdir -p $REFERENCE_DATA_DIR
zenodo_get  $REFERENCE_DATA_ID -o $REFERENCE_DATA_DIR --continue-on-error --retry 2 --pause 30