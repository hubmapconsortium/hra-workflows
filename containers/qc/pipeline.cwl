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
  - --output-matrix
  - filtered_matrix.h5ad

inputs:
  matrix:
    type: File
    label: Data to run QC on in h5ad format
    inputBinding:
      position: 0
  organ:
    type: string
    label: Organ uberon id in format 'UBERON:1234'
    inputBinding:
      prefix: --organ
  options: ./options.yml#options

outputs:
  filtered_matrix:
    type: File
    outputBinding:
      glob: filtered_matrix.h5ad
    doc: "Matrix produced by the QC step for downstream annotation."
  report:
    type: File
    outputBinding:
      glob: report.json
  qc_results:
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

  filtered_qc_matrix:
    type: File?
    outputBinding:
      glob: "qc_results/filtered_data.h5ad"
    doc: "Filtered AnnData file containing only high-quality cells (if --filter)."
