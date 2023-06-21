#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0

requirements:
  DockerRequirement:
    dockerPull: ghcr.io/hubmapconsortium/hra-workflows/celltypist:main

inputs:
  matrix:
    type: File
    doc: Data to annotate
    inputBinding:
      position: 0
  referenceData:
    type: File?
    doc: Not used by celltypist
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
        fields:
          existingAnnotationsColumn:
            type: string?
            doc: Column with existing annotation to compare predictions against
            inputBinding:
              prefix: --existing-annotations-column=
              separate: false

outputs:
  annotations:
    type: File
    outputBinding:
      glob: annotations.csv
  annotated_matrix:
    type: File
    outputBinding:
      glob: annotated_matrix.h5ad
