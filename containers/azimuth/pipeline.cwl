#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0

requirements:
  DockerRequirement:
    dockerPull: ghcr.io/hubmapconsortium/hra-workflows/azimuth:main

inputs:
  matrix:
    type: File
    doc: Data to annotate
    inputBinding:
      position: 0
  referenceData:
    type: File
    doc: Reference data set
    inputBinding:
      prefix: --reference-data=
      separate: false
  organ:
    type: string
    doc: Organ uberon id in format 'UBERON:1234'
    inputBinding:
      prefix: --organ=
      separate: false
  options:
    type:
      - "null"
      - type: record
        fields: {}

outputs:
  annotations:
    type: File
    outputBinding:
      glob: annotations.csv
  annotated_matrix:
    type: File
    outputBinding:
      glob: annotated_matrix.h5ad
