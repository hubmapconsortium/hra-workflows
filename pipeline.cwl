#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.0

requirements:
  MultipleInputFeatureRequirement: {}
  SubworkflowFeatureRequirement: {}

inputs:
  matrix:
    type: File
  organ:
    type: string

outputs:
  output_dirs:
    type: Directory[]
    outputSource: collect/dirs

steps:
  azimuth:
    run: containers/azimuth/pipeline.cwl
    in:
      matrix: matrix
      organ: organ
    out: [annotations, annotated_matrix]

  celltypist:
    run: containers/celltypist/pipeline.cwl
    in:
      matrix: matrix
      organ: organ
    out: [annotations, annotated_matrix]

  popv:
    run: containers/popv/pipeline.cwl
    in:
      matrix: matrix
      organ: organ
    out: [annotations, annotated_matrix]

  collect:
    run:
      doc: Collect output files from each algorithm into separate output directories
      class: ExpressionTool
      requirements:
        - class: InlineJavascriptRequirement
          expressionLib:
            - |
              function inputsToDirs(inputs) {
                var dirs = [];
                for (var group in inputs) {
                  dirs.push({
                    class: "Directory",
                    basename: group,
                    listing: inputs[group]
                  });
                }

                return dirs;
              }

      inputs:
        azimuth: File[]
        celltypist: File[]
        popv: File[]

      outputs:
        dirs: Directory[]

      expression: |
        $({ dirs: inputsToDirs(inputs) })

    in:
      azimuth: [azimuth/annotations, azimuth/annotated_matrix]
      celltypist: [celltypist/annotations, celltypist/annotated_matrix]
      popv: [popv/annotations, popv/annotated_matrix]
    out: [dirs]
