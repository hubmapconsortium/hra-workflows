#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0

requirements:
  DockerRequirement:
    dockerPull: ghcr.io/hubmapconsortium/hra-workflows/extract-summary:main

inputs:
  annotations:
    type: File
    doc: Annotations csv file
    inputBinding:
      position: 0
  annotationMethod:
    type: string
    doc: Method used to extract annotations
    inputBinding:
      prefix: --annotation-method=
      separate: false
  cellLabelColumn:
    type: string
    doc: Cell label column. Used for grouping if --cell-id-column is not provided.
    inputBinding:
      prefix: --cell-label-column=
      separate: false
  cellIdColumn:
    type: string?
    doc: Optional cell id column. Groups by label if not provided.
    inputBinding:
      prefix: --cell-id-column=
      separate: false
  cellSource:
    type: string?
    doc: Cell source. Must be an IRI.
    inputBinding:
      prefix: --cell-source=
      separate: false
  jsonldContext:
    type: File?
    doc: Base jsonld context to add summary to.
    inputBinding:
      prefix: --jsonld-context=
      separate: false

outputs:
  summary:
    type: File
    outputBinding:
      glob: summary.jsonld
