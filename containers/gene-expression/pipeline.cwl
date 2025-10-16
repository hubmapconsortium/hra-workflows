#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0

requirements:
  DockerRequirement:
    dockerPull: ghcr.io/hubmapconsortium/hra-workflows/gene-expression:main
  SchemaDefRequirement:
    types:
      - $import: ./options.yml

baseCommand: python
arguments:
  - /main.py

inputs:
  matrix:
    type: File
    label: Data to get gene expression for in h5ad format
    inputBinding:
      position: 0
  options: ./options.yml#options

outputs:
  matrix_with_gene_expr:
    type: File?
    outputBinding:
      glob: matrix_with_gene_expr.h5ad
  report:
    type: File?
    outputBinding:
      glob: report.json
