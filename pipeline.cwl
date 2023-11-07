#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2

requirements:
  MultipleInputFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  SubworkflowFeatureRequirement: {}
  SchemaDefRequirement:
    types:
      - $import: ./containers/azimuth/options.yml
      - $import: ./containers/celltypist/options.yml
      - $import: ./containers/popv/options.yml
      - $import: ./containers/crosswalking/options.yml
      - $import: ./containers/gene-expression/options.yml
      - $import: ./containers/extract-summary/options.yml

inputs:
  matrix: File
  organ: string
  algorithms:
    type:
      type: array
      items:
        type: record
        fields:
          azimuth: ./containers/azimuth/options.yml#options?
          celltypist: ./containers/celltypist/options.yml#options?
          popv: ./containers/popv/options.yml#options?
          crosswalk: ./containers/crosswalking/options.yml#options?
          geneExpression: ./containers/gene-expression/options.yml#options?
          summarize: ./containers/extract-summary/options.yml#options?
          directory: string?

outputs:
  directories:
    type: Directory[]
    outputSource: runEach/directory

steps:
  runEach:
    run: ./steps/run-one.cwl
    scatter: algorithm
    in:
      matrix: matrix
      organ: organ
      algorithm: algorithms
    out: [directory]
