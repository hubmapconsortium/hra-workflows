library(Azimuth)
library(Seurat)
library(SeuratData)
library(SeuratDisk)

args <- commandArgs(trailingOnly = TRUE)

matrix_path <- args[1]
reference <- args[2]

# Annotate
output_data <- RunAzimuth(matrix_path, reference=reference)

# Save and convert to h5ad
SaveH5Seurat(output_data, 'result.h5seurat')
Convert('result.h5seurat', dest='h5ad')
