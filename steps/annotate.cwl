#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2

requirements:
  InlineJavascriptRequirement: {}
  MultipleInputFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}
  SchemaDefRequirement:
    types:
      - $import: ../containers/azimuth/options.yml
      - $import: ../containers/celltypist/options.yml
      - $import: ../containers/popv/options.yml
      - $import: ../containers/pan-human-azimuth/options.yml
      - $import: ../containers/frmatch/options.yml

inputs:
  matrix: File
  organ: string
  algorithm:
    type:
      type: record
      fields:
          azimuth: ../containers/azimuth/options.yml#options?
          celltypist: ../containers/celltypist/options.yml#options?
          popv: ../containers/popv/options.yml#options?
          pan-human-azimuth: ../containers/pan-human-azimuth/options.yml#options?
          frmatch: ../containers/frmatch/options.yml#options?

outputs:
  annotations:
    type: File
    outputSource:
      - azimuth/annotations
      - celltypist/annotations
      - popv/annotations
      - pan-human-azimuth/annotations
      - frmatch/annotations
    pickValue: first_non_null
  annotated_matrix:
    type: File
    outputSource:
      - azimuth/annotated_matrix
      - celltypist/annotated_matrix
      - popv/annotated_matrix
      - pan-human-azimuth/annotated_matrix
      - frmatch/annotated_matrix
    pickValue: first_non_null
  report:
    type: File
    outputSource:
      - azimuth/report
      - celltypist/report
      - popv/report
      - pan-human-azimuth/report
      - frmatch/report
    pickValue: first_non_null

steps:
  azimuth:
    run: ../containers/azimuth/pipeline.cwl
    when: $(!!inputs.options)
    in:
      matrix: matrix
      organ: organ
      options:
        source: algorithm
        valueFrom: $(inputs.options.azimuth || null)
    out: [annotations, annotated_matrix, report]
  
  celltypist:
    run: ../containers/celltypist/pipeline.cwl
    when: $(!!inputs.options)
    in:
      matrix: matrix
      organ: organ
      options:
        source: algorithm
        valueFrom: $(inputs.options.celltypist || null)
    out: [annotations, annotated_matrix, report]

  popv:
    run: ../containers/popv/pipeline.cwl
    when: $(!!inputs.options)
    in:
      matrix: matrix
      organ: organ
      options:
        source: algorithm
        valueFrom: $(inputs.options.popv || null)
    out: [annotations, annotated_matrix, report]

  pan-human-azimuth:
    run: ../containers/pan-human-azimuth/pipeline.cwl
    when: $(!!inputs.options)
    in:
      matrix: matrix
      organ: organ
      options:
        source: algorithm
        valueFrom: $(inputs.options['pan-human-azimuth'] || null)
    out: [annotations, annotated_matrix, report]
  
  frmatch:
    run: ../containers/frmatch/pipeline.cwl
    when: $(!!inputs.options)
    in:
      matrix: matrix
      organ: organ
      options:
        source: algorithm
        valueFrom: $(inputs.options.frmatch || null)
    out: [annotations, annotated_matrix, report]
