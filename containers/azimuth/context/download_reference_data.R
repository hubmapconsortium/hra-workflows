library(rjson)
library(Seurat)
library(SeuratData)

args <- commandArgs(trailingOnly = TRUE)
organ_mapping_file <- args[1]
output_dir <- args[2]

mapping <- fromJSON(file = organ_mapping_file)
references <- unlist(unique(mapping))
InstallData(references, destdir=output_dir)
