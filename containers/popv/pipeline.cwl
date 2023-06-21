#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0

requirements:
  DockerRequirement:
    dockerPull: ghcr.io/hubmapconsortium/hra-workflows/popv:main

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
        fields:
          cellOntologyDirectory:
            type: string?
            inputBinding:
              prefix: --cell-ontology-dir=
              separate: false
          queryLabelsKey:
            type: string?
            inputBinding:
              prefix: --query-labels-key=
              separate: false
          queryBatchKey:
            type: string?
            inputBinding:
              prefix: --query-batch-key=
              separate: false
          referenceLabelsKey:
            type: string?
            inputBinding:
              prefix: --ref-labels-key=
              separate: false
          referenceBatchKey:
            type: string?
            inputBinding:
              prefix: --ref-batch-key=
              separate: false
          unknownLabelsKey:
            type: string?
            inputBinding:
              prefix: --unknown-labels-key=
              separate: false
          samplesPerLabel:
            type: int?
            inputBinding:
              prefix: --samples-per-label=
              separate: false
          sampleSize:
            type: int?
            inputBinding:
              prefix: --sample-size=
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
