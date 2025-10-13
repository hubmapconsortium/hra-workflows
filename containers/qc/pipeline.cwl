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
  - "--log-level"
  - "40"

inputs:
  matrix:
    type: File
    label: Data to run QC on in h5ad format
    inputBinding:
      position: 0
  options:
    type: "./options.yml#options?"

outputs:
  outputDir:
    type: Directory
    outputBinding:
      glob: 'qc_results'
    doc: "QC results"

  qcPerCellCsv:
    type: File
    outputBinding:
      glob: "qc_results/qc_per_cell.csv"
    doc: "Per-cell QC metrics table."

  qcSummaryCsv:
    type: File
    outputBinding:
      glob: "qc_results/qc_summary.csv"
    doc: "Summary QC statistics across cells."

  qcSummaryJson:
    type: File
    outputBinding:
      glob: "qc_results/qc_summary.json"
    doc: "JSON summary of QC results and thresholds."

  filteredH5ad:
    type: File?
    outputBinding:
      glob: "qc_results/filtered_data.h5ad"
    doc: "Filtered AnnData file containing only high-quality cells (if --filter)."
