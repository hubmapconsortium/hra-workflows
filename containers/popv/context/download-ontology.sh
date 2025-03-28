#!/bin/bash

ONTOLOGY_DIR=./ontology
ONTOLOGY_DIR_URL=https://raw.githubusercontent.com/YosefLab/PopV/2d29c9a290d2015ec65ef0ef9f0e6b6d2277e7bb/resources/ontology

mkdir -p $ONTOLOGY_DIR
cd $ONTOLOGY_DIR
rm -f *

wget $ONTOLOGY_DIR_URL/cl.obo
wget $ONTOLOGY_DIR_URL/cl.ontology
wget $ONTOLOGY_DIR_URL/cl.ontology.nlp.emb
