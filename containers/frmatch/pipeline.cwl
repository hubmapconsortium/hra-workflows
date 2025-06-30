#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0

requirements:
  DockerRequirement:
    dockerPull: ghcr.io/hubmapconsortium/hra-workflows/frmatch:main
  SchemaDefRequirement:
    types:
      - $import: ./options.yml

baseCommand: python
arguments:
  - /main.py

inputs:
  matrix:
    type: File
    label: Data to annotate in h5ad format
    inputBinding:
      position: 0
  organ:
    type: string
    label: Organ uberon id in format 'UBERON:1234'
    inputBinding:
      prefix: --organ
  options: ./options.yml#options

outputs:
  annotations:
    type: File
    outputBinding:
      glob: annotations.csv
  annotated_matrix:
    type: File
    outputBinding:
      glob: annotated_matrix.h5ad
  report:
    type: File
    outputBinding:
      glob: report.json
