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
    doc: Annotated matrix h5ad
    inputBinding:
      position: 0
  options: ./options.yml#options

outputs:
  annotations:
    type: File
    outputBinding:
      glob: annotations.csv.gz
  summary:
    type: File
    outputBinding:
      glob: summary.jsonld
