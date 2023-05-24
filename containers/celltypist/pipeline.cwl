#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0

requirements:
  DockerRequirement:
    dockerPull: ghcr.io/hubmapconsortium/hra-workflows/celltypist:main

inputs:
  matrix:
    type: File
    inputBinding:
      position: 1
  model:
    type: string
    inputBinding:
      position: 2
      prefix: --model=
      separate: false
  existingAnnotationsColumn:
    type: string?
    inputBinding:
      position: 3
      prefix: --existing-annotations-column=
      separate: false

outputs:
  annotations:
    type: File
    outputBinding:
      glob: annotations.csv
