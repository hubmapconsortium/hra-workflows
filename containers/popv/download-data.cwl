#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2

requirements:
  DockerRequirement:
    dockerPull: ghcr.io/hubmapconsortium/hra-workflows/popv:main
  NetworkAccess:
    networkAccess: true

baseCommand: /bin/bash
arguments:
  - /download-data.sh

inputs:
  modelsZenodoId:
    type: string
    label: Zenodo id for models collection
    default: '7580707'
    inputBinding:
      position: 0
  referenceDataZenodoId:
    type: string
    label: Zenodo id for reference data collection
    default: '7587774'
    inputBinding:
      position: 1
  outputDirectory:
    type: string
    label: Output directory for models and reference data
    default: ./popv
    inputBinding:
      position: 2

outputs:
  data:
    type: Directory
    outputBinding:
      glob: $(inputs.outputDirectory)
      loadListing: deep_listing
