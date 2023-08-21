#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2

requirements:
  MultipleInputFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  SubworkflowFeatureRequirement: {}
  SchemaDefRequirement:
    types:
      - $import: ./containers/preprocess/options.yml
      - $import: ./containers/azimuth/options.yml
      - $import: ./containers/celltypist/options.yml
      - $import: ./containers/popv/options.yml
      - $import: ./containers/extract-summary/options.yml

inputs:
  matrix: File
  organ: string
  preprocessing: ./containers/preprocess/options.yml#options
  algorithms:
    type:
      type: array
      items:
        type: record
        fields:
          azimuth: ./containers/azimuth/options.yml#options?
          celltypist: ./containers/celltypist/options.yml#options?
          popv: ./containers/popv/options.yml#options?
          extract: ./containers/extract-summary/options.yml#options
          directory: string?

outputs:
  directories:
    type: Directory[]
    outputSource: runEach/directory

steps:
  preprocess:
    run: ./containers/preprocess/pipeline.cwl
    in:
      matrix: matrix
      options: preprocessing
    out: [preprocessed_matrix]
  runEach:
    run: ./steps/run-one.cwl
    scatter: algorithm
    in:
      matrix: preprocess/preprocessed_matrix
      organ: organ
      algorithm: algorithms
    out: [directory]
