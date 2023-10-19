library(Azimuth)
library(Seurat)
library(SeuratData)
library(SeuratDisk)

args <- commandArgs(trailingOnly = TRUE)

organ <- args[1]
matrix_path <- args[2]
# dest_matrix_path = gsub("data","annotated-data",matrix_path)

# data <- LoadFileInput(path=matrix_path)
print('start')
output_data <- RunAzimuth(matrix_path, reference = organ)
print('done')
SaveH5Seurat(output_data, 'result.h5seurat')
Convert('result.h5seurat', dest='h5ad')
