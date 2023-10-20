library(rjson)
library(Seurat)
library(SeuratData)

args <- commandArgs(trailingOnly = TRUE)
organ_mapping_file <- args[1]
output_dir <- args[2]

# Load unique reference organs
mapping <- fromJSON(file = organ_mapping_file)
references <- unlist(unique(mapping))

# Download and install data
options(timeout=60 * 60) # Probably overkill but the default of 60s is to low for some of the datasets
InstallData(references, lib=output_dir)
