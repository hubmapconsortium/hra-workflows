#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2

requirements:
  MultipleInputFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}
  InlineJavascriptRequirement:
    expressionLib:
      - $include: ./js/options-util.js

  SchemaDefRequirement:
    types:
      - $import: ../containers/azimuth/options.yml
      - $import: ../containers/celltypist/options.yml
      - $import: ../containers/popv/options.yml
      - $import: ../containers/crosswalking/options.yml
      - $import: ../containers/gene-expression/options.yml
      - $import: ../containers/extract-summary/options.yml

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
          crosswalk: ../containers/crosswalking/options.yml#options?
          geneExpression: ../containers/gene-expression/options.yml#options?
          summarize: ../containers/extract-summary/options.yml#options?
          directory: string?

outputs:
  directory:
    type: Directory
    outputSource: collect/directory

steps:
  annotate:
    run: ./annotate.cwl
    in:
      matrix: matrix
      organ: organ
      algorithm: algorithm
    out: [annotations, annotated_matrix, report]

  check_result:
    run: ./check_annotation_report.cwl
    in:
      report: annotate/report
      matrix: annotate/annotated_matrix
    out: [matrix_or_null]

  crosswalk:
    run: ../containers/crosswalking/pipeline.cwl
    when: $(!!inputs.matrix)
    in:
      matrix: check_result/matrix_or_null
      options:
        source: algorithm
        valueFrom: $(self.crosswalk || {})
    out: [matrix_with_crosswalking]

  gene_expression:
    run: ../containers/gene-expression/pipeline.cwl
    when: $(!!inputs.matrix)
    in:
      matrix: crosswalk/matrix_with_crosswalking
      options:
        source: algorithm
        valueFrom: $(self.geneExpression || {})
    out: [matrix_with_gene_expr]
  
  summarize:
    run: ../containers/extract-summary/pipeline.cwl
    when: $(!!inputs.matrix)
    in:
      matrix: gene_expression/matrix_with_gene_expr
      options:
        source: algorithm
        valueFrom: $(self.summarize || {})
    out: [summary, annotations]

  collect:
    run: ./collect-files.cwl
    in:
      files:
        source: [summarize/summary, summarize/annotations, annotate/report]
        pickValue: all_non_null
      outputDirectory:
        source: algorithm
        valueFrom: $(self.directory || {})
    out: [directory]
