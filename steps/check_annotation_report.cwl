#!/usr/bin/env cwl-runner
class: ExpressionTool
cwlVersion: v1.2

requirements:
  InlineJavascriptRequirement: {}

inputs:
  report:
    type: File
    loadContents: true
  matrix: File

outputs:
  matrix_or_null: File?

expression: |
  ${
    var isSuccess = /"status": "success"/g.test(inputs.report.contents);
    return { matrix_or_null: isSuccess ? inputs.matrix : null };
  }
