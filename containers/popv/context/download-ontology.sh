#!/bin/bash

ONTOLOGY_DIR=./ontology
ONTOLOGY_DIR_URL=https://raw.githubusercontent.com/YosefLab/PopV/main/ontology

mkdir -p $ONTOLOGY_DIR
cd $ONTOLOGY_DIR
rm -f *

wget $ONTOLOGY_DIR_URL/cl.obo
wget $ONTOLOGY_DIR_URL/cl.ontology
wget $ONTOLOGY_DIR_URL/cl.ontology.nlp.emb
