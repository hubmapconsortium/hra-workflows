#!/bin/bash
shopt -s extglob
set -e

ALL_CONTAINERS=$(basename -a ./containers/*)
CONTAINERS=${@:-$ALL_CONTAINERS}

for NAME in $CONTAINERS; do
  cp -r ./src/ ./containers/$NAME/context/
  docker build -t ghcr.io/hubmapconsortium/hra-workflows/$NAME:main ./containers/$NAME/
done
