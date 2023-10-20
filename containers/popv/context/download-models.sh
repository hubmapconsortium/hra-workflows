#!/bin/bash
set -e

MODELS_ID=${1:?"A zenodo models id must be provided to download!"}
MODELS_DIR=${2:-"./popv/models"}

mkdir -p $MODELS_DIR
zenodo_get $MODELS_ID -o $MODELS_DIR

for ARCHIVE in $MODELS_DIR/*.tar.gz; do
  MODEL=$(basename -s .tar.gz $ARCHIVE)
  DIR="$MODELS_DIR/$MODEL"

  if [[ ! -d $DIR || -z $(ls -A $DIR) ]]; then
    mkdir -p $DIR
    tar zx -C $DIR -f $ARCHIVE
  fi
done
