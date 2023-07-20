#!/usr/bin/env cwl-runner
class: ExpressionTool
cwlVersion: v1.2

requirements:
  InlineJavascriptRequirement: {}
  MultipleInputFeatureRequirement: {}

inputs:
  files: File[]
  outputDirectory: string
outputs:
  directory: Directory

expression: |
  ${
    var parts = inputs.outputDirectory.split('/');
    var listing = inputs.files;
    for (var i = parts.length - 1; i >= 0; --i) {
      listing = [{
        class: 'Directory',
        basename: parts[i],
        listing: listing
      }];
    }

    return { directory: listing[0] }
  }
