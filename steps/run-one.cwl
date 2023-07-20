#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2

requirements:
  MultipleInputFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}
  InlineJavascriptRequirement:
    expressionLib:
      - |
        function getDirectoryName(algorithm) {
          var nonAlgorithms = ['extract', 'directory'];
          var algorithms = Object.keys(algorithm).filter(function(key) {
            return !nonAlgorithms.includes(key) && !!algorithm[key];
          });
          var candidates = [].concat(algorithm.directory, algorithms, '.');
          return candidates.find(function(value) { return !!value });
        }

  SchemaDefRequirement:
    types:
      - $import: ../containers/azimuth/options.yml
      - $import: ../containers/celltypist/options.yml
      - $import: ../containers/popv/options.yml
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
          extract: ../containers/extract-summary/options.yml#options
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
    out: [annotations, annotated_matrix]
  
  extract:
    run: ../containers/extract-summary/pipeline.cwl
    in:
      annotations: annotate/annotations
      options:
        source: algorithm
        valueFrom: $(self.extract)
    out: [summary]

  collect:
    run: ./collect-files.cwl
    in:
      files: [annotate/annotations, annotate/annotated_matrix, extract/summary]
      outputDirectory:
        source: algorithm
        valueFrom: $(getDirectoryName(self))
    out: [directory]