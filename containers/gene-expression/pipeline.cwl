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
  - --output-report
  - report_gene_expr.json
  - --output-gene-expr-json
  - gene_expr.json

inputs:
  matrix:
    type: File
    label: Data to get gene expression for in h5ad format
    inputBinding:
      position: 0
  options: ./options.yml#options?

outputs:
  gene_expr_json:
    type: File
    outputBinding:
      glob: "gene_expr.json"
  report:
    type: File?
    outputBinding:
      glob: "report_gene_expr.json"
  