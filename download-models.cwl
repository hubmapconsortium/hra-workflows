#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2

requirements:
  InlineJavascriptRequirement: {}
  MultipleInputFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}

inputs:
  outputDirectory:
    type: string
    label: Output directory for model data
    default: ./models

outputs:
  data:
    type: Directory
    outputSource: collectFiles/directory

steps:
  downloadAzimuth:
    run: ./containers/azimuth/download-data.cwl
    in: []
    out: [data]

  downloadPopv:
    run: ./containers/popv/download-data.cwl
    in: []
    out: [data]

  collectFiles:
    run: ./steps/collect-files.cwl
    in:
      files: [downloadAzimuth/data, downloadPopv/data]
      outputDirectory: outputDirectory
    out: [directory]
