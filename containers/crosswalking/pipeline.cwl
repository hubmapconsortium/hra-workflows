#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0

requirements:
  DockerRequirement:
    dockerPull: ghcr.io/hubmapconsortium/hra-workflows/crosswalking:main
  SchemaDefRequirement:
    types:
      - $import: ./options.yml

baseCommand: python
arguments:
  - /main.py

inputs:
  matrix:
    type: File
    label: Data to get crosswalking for in h5ad format
    inputBinding:
      position: 0
  options: ./options.yml#options

outputs:
  annotations:
    type: File
    outputBinding:
      glob: annotations.csv.gz
  matrix_with_crosswalking:
    type: File
    outputBinding:
      glob: matrix_with_crosswalking.h5ad
