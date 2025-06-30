#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2

requirements:
  DockerRequirement:
    dockerPull: ghcr.io/hubmapconsortium/hra-workflows/frmatch:main
  NetworkAccess:
    networkAccess: true

baseCommand: /bin/bash
arguments:
  - /download-data.sh

inputs:
  referenceDataZenodoId:
    type: string
    label: Zenodo id for reference data collection
    default: '15750919'
    inputBinding:
      position: 0
  outputDirectory:
    type: string
    label: Output directory for models and reference data
    default: ./frmatch
    inputBinding:
      position: 1

outputs:
  data:
    type: Directory
    outputBinding:
      glob: $(inputs.outputDirectory)
      loadListing: deep_listing
