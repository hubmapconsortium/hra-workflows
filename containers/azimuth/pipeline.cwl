#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.0

requirements:
  SubworkflowFeatureRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: secondary_analysis.h5ad
        entry: ''
      - entryname: version_metadata.json
        entry: ''
      - entryname: annotations.csv
        entry: ''

inputs:
  matrix:
    type: File
  # Generalize other inputs
  reference:
    # TODO figure out how to make an enum
    # Supported values: LK, RK, LL, RL, HT
    type: string
  secondary_analysis_matrix:
    # TODO figure out what this does. Currently I use the same file as matrix for this input
    type: File

outputs:
  # Output other files? (version metadata and annotated secondary analysis)
  annotations:
    type: File
    outputSource: azimuth/annotations_csv

steps:
  azimuth:
    run: https://raw.githubusercontent.com/hubmapconsortium/azimuth-annotate/main/steps/azimuth-annotate.cwl
    in:
      matrix:
        source: matrix
      reference:
        source: reference
      secondary_analysis_matrix:
        source: secondary_analysis_matrix
    out: [annotations_csv]
