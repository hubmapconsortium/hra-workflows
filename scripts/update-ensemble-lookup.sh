#!/bin/bash

URL=https://apps.humanatlas.io/api/grlc/hra/ensembl-lookup.csv

curl -o src/assets/ensemble-lookup.csv $URL
