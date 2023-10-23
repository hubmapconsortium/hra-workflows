#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2

requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}

inputs:
  outputDirectory:
    type: string
    label: Output directory for model data
    default: ./models

outputs:
  data:
    type: Directory[]
    outputSource: [downloadAzimuth/data, downloadPopv/data]

steps:
  downloadAzimuth:
    run: ./containers/azimuth/download-data.cwl
    in:
      outputDirectory:
        source: outputDirectory
        valueFrom: $(self + '/azimuth')
    out: [data]

  downloadPopv:
    run: ./containers/popv/download-data.cwl
    in:
      outputDirectory:
        source: outputDirectory
        valueFrom: $(self + '/popv')
    out: [data]
