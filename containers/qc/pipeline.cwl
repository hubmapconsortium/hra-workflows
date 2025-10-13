#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0

requirements:
  DockerRequirement:
    dockerPull: ghcr.io/hubmapconsortium/hra-workflows/qc:main
  SchemaDefRequirement:
    types:
      - $import: ./options.yml

baseCommand: python
arguments:
  - /main.py
  - --output
  - "$(inputs.options.outputDir)"
  - "--log-level"
  - "40"

inputs:
  matrix:
    type: File
    label: Data to annotate in h5ad format
    inputBinding:
      position: 0
  options: ./options.yml#options

outputs:
  qcPerCellCsv:
    type: File
    outputBinding:
      glob: "$(inputs.outputDir)/qc_per_cell.csv"
    doc: "Per-cell QC metrics table."

  qcSummaryCsv:
    type: File
    outputBinding:
      glob: "$(inputs.outputDir)/qc_summary.csv"
    doc: "Summary QC statistics across cells."

  qcSummaryJson:
    type: File
    outputBinding:
      glob: "$(inputs.outputDir)/qc_summary.json"
    doc: "JSON summary of QC results and thresholds."

  filteredH5ad:
    type: File?
    outputBinding:
      glob: "$(inputs.outputDir)/filtered_data.h5ad"
    doc: "Filtered AnnData file containing only high-quality cells (if --filter)."
