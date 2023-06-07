#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0

requirements:
  DockerRequirement:
    dockerPull: ghcr.io/hubmapconsortium/hra-workflows/popv:main

inputs:
  matrix:
    type: File
    inputBinding:
      position: 1
  organ:
    type: string
    inputBinding:
      position: 2
      prefix: --organ=
      separate: false
  referenceData:
    type: File?
    inputBinding:
      position: 3
      prefix: --reference-organ=
      separate: false
  cellOntologyDirectory:
    type: string?
    inputBinding:
      position: 4
      prefix: --cell-ontology-dir=
      separate: false
  queryLabelsKey:
    type: string?
    inputBinding:
      position: 5
      prefix: --query-labels-key=
      separate: false
  queryBatchKey:
    type: string?
    inputBinding:
      position: 6
      prefix: --query-batch-key=
      separate: false
  referenceLabelsKey:
    type: string?
    inputBinding:
      position: 7
      prefix: --ref-labels-key=
      separate: false
  referenceBatchKey:
    type: string?
    inputBinding:
      position: 8
      prefix: --ref-batch-key=
      separate: false
  unknownLabelsKey:
    type: string?
    inputBinding:
      position: 9
      prefix: --unknown-labels-key=
      separate: false
  samplesPerLabel:
    type: int?
    inputBinding:
      position: 10
      prefix: --samples-per-label=
      separate: false
  sampleSize:
    type: int?
    inputBinding:
      position: 11
      prefix: --sample-size=
      separate: false

outputs:
  annotations:
    type: File
    outputBinding:
      glob: annotations.csv
