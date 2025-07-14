# hra-workflows
HRA workflows for annotating h5ad using different tools.

## Table of Content
- [Requirements](#requirements)
- [Building docker images](#building-the-docker-images-locally)
- [Running the annotation tools](#running-the-annotation-tools)
  - [Downloading models](#1-download-annotation-tool-models)
  - [Creating jobs](#2-create-a-job-file)
  - [Running jobs](#3-running-the-job)
- [Adding custom annotation tools](#adding-a-new-annotation-tool)
  - [Folder structure](#1-tool-folder-structure)
  - [Updating top level pipeline](#2-updating-top-level-steps)

#
## Requirements
- [cwl-runner](https://github.com/common-workflow-language/cwltool)
- Docker

## Building the docker images locally
Docker images can be build locally by running `./scripts/build-containers.sh`. By default the script will build all containers when run. To build individual images or a set of images provide the container names as arguments to the build script, ex. `./scripts/build-containers.sh azimuth gene-expression`.

## Running the annotation tools

### 1. Download annotation tool models
Download model data by running `cwl-runner download-models.cwl`.

### 2. Create a job file
The first step is to create a job file that will specify inputs to the pipeline. The file can be written as either a json or yaml file.

#### Example job.yml running Azimuth
```yaml
matrix:
  class: File
  path: path/to/data.h5ad
organ: UBERON:0002048 # Uberon id for lung
algorithms:
  # Algorithm specific options are documented in the container's options.yml
  - azimuth:
      referenceDataDir:
        class: Directory
        path: path/to/models/directory
```

### 3. Running the job
After creating a job file running the annotation tools is as simple as running `cwl-runner pipeline.cwl my-job.yml` (replace `my-job.yml` with your job file).

## Adding a new annotation tool

### 1. Tool folder structure
An annotation tool generally has the following file structure:
```
containers/
  my-annotation-tool/
    Dockerfile
    options.yml
    pipeline.cwl
    download-data.cwl (optional)
    context/*
      code and assets...
```

Where each file should perform the following function:
- `Dockerfile`
  - Instructions for building a docker image.
- `options.yml`
    - Cwl definition of tool specific options.
- `pipeline.cwl`
  - Main cwl pipeline for running the tool.
  - 3 inputs: "matrix", "organ", and "options"
  - 3 outputs:  "annotations", "annotated_matrix", and "report".
- `download-data.cwl` (optional)
  - Download models and other data required for running the tool.
  - Implement this pipeline when the model data is to large to embed directly in the docker image.
- `context/*`
  - Directory containing the code and assets implementing the tool.

### 2. Updating top level steps
After implementing a new algorithm a few changes have to be made to enable the tool from the main pipeline. The files that have to be updated are: `pipeline.cwl`, `./steps/annotate.cwl`, and `./steps/run-one.cwl`. After adding the new tool to the top level pipeline it can be used by specifying the tool in a job file.

#### Example job for a new annotation tool
```yaml
matrix:
  class: File
  path: path/to/data.h5ad
organ: UBERON:0002048 # Uberon id for lung
algorithms:
  - my-annotation-tool:
      # Options specific to my-annotation-tool
      option1: value1
      option2: value2
      ...
```
