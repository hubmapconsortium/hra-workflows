#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0

requirements:
  DockerRequirement:
    dockerPull: ghcr.io/hubmapconsortium/hra-workflows/preprocess:main
  SchemaDefRequirement:
    types:
      - $import: ./options.yml

baseCommand: python /main.py

inputs:
  matrix:
    type: File
    label: Data to preprocess in h5ad format
    inputBinding:
      position: 0
  options: ./options.yml#options

outputs:
  preprocessed_matrix:
    type: File
    outputBinding:
      glob: preprocessed_matrix.h5ad
