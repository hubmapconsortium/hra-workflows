#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0

requirements:
  DockerRequirement:
    dockerPull: ghcr.io/hubmapconsortium/hra-workflows/nsforest:main
  SchemaDefRequirement:
    types:
      - $import: ./options.yml

baseCommand: python
arguments:
  - /main.py
  - --output-report
  - report_nsforest.json

inputs:
  matrix:
    type: File
    label: Data to get NSForest markers for in h5ad format
    inputBinding:
      position: 0
  options: ./options.yml#options

outputs:
  nsforest_gene_expr_json:
    type: File
    outputBinding:
      glob: nsforest_gene_expr.json
  report:
    type: File?
    outputBinding:
      glob: report_nsforest.json
