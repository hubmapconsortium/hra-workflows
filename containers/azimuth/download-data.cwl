#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2

requirements:
  DockerRequirement:
    dockerPull: ghcr.io/hubmapconsortium/hra-workflows/azimuth:main
  NetworkAccess:
    networkAccess: true

baseCommand: /bin/bash
arguments:
  - /download-data.sh

inputs:
  outputDirectory:
    type: string
    label: Output directory for reference data
    default: ./azimuth
    inputBinding:
      position: 0
  organMappingFile:
    type: File?
    label: Organ mapping json file
    inputBinding:
      position: 1

outputs:
  data:
    type: Directory
    outputBinding:
      glob: $(inputs.outputDirectory)
      loadListing: deep_listing
