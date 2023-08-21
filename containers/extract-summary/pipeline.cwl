#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0

requirements:
  DockerRequirement:
    dockerPull: ghcr.io/hubmapconsortium/hra-workflows/extract-summary:main
  SchemaDefRequirement:
    types:
      - $import: ./options.yml

baseCommand: python /main.py

inputs:
  annotations:
    type: File
    doc: Annotations csv file
    inputBinding:
      position: 0
  options: ./options.yml#options

outputs:
  summary:
    type: File
    outputBinding:
      glob: summary.jsonld
