#!/bin/bash

# Defined by user: MODELS_DIR and ZENODO_MODELS_ID
MODELS_TMP_DIR=$MODELS_DIR/tmp

mkdir -p $MODELS_DIR $MODELS_TMP_DIR

pip install 'zenodo_get@git+https://github.com/dvolgyes/zenodo_get'
zenodo_get $ZENODO_MODELS_ID -o $MODELS_TMP_DIR

for archive in $MODELS_TMP_DIR/*.tar.gz; do
  NAME=`basename $archive .tar.gz`
  mkdir -p $MODELS_DIR/${NAME}
  tar zx -C $MODELS_DIR/${NAME} -f $archive
done

rm -rf $MODELS_TMP_DIR
