#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0

requirements:
  DockerRequirement:
    dockerPull: ghcr.io/hubmapconsortium/hra-workflows/extract-summary:main
  SchemaDefRequirement:
    types:
      - $import: ./options.yml

baseCommand: python
arguments:
  - /main.py

inputs:
  matrix:
    type: File
    doc: Matrix h5ad for cell labels and counts
    inputBinding:
      position: 0
  gene_expr_json:
    type: File
    doc: Gene expression data from JSON file
    inputBinding:
      prefix: --gene-expr-json
  nsforest_gene_expr_json:
    type: File
    doc: NSForest gene expression data from JSON file
    inputBinding:
      prefix: --nsforest-gene-expr-json
  options: ./options.yml#options

outputs:
  summary:
    type: File
    outputBinding:
      glob: summary.jsonld
