#!/usr/bin/env bash

### In linux terminal or web browser, download published canine data -- put in "external_data" dir
# wget https://zenodo.org/record/11434076/files/AllCells_duod_annotated.rds
# wget https://zenodo.org/record/10891255/files/canine_naive_n6_annotated.rds
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225599/suppl/GSE225599_final_dataSet_HvO.rds.gz
# gunzip GSE225599_final_dataSet_HvO.rds.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225599/suppl/GSE225599_final_dataSet_H.rds.gz
# gunzip GSE225599_final_dataSet_H.rds.gz

## Load libraries
library(Seurat)
library(SeuratDisk)
## Prep canine naive OS reference
seu.obj <- readRDS("../external_data/canine_naive_n6_annotated.rds")
seu.obj$celltype.l3 <- as.character(seu.obj$celltype.l3)
SaveH5Seurat(seu.obj, filename = "../external_data/canine_naive_n6_annotated.h5Seurat", verbose = TRUE, overwrite = T)
Convert("../external_data/canine_naive_n6_annotated.h5Seurat", dest = "h5ad", overwrite = T)
## Prep canine OS query dataset
seu.obj <- readRDS("../internal_data/20230504_naiveTX_n6_n6_QCfiltered_1_2500_res0.8_dims45_dist0.3_neigh30_S3.rds")
SaveH5Seurat(seu.obj, filename = "../internal_data/canine_OS_query.h5Seurat", verbose = TRUE, overwrite = T)
Convert("../internal_data/canine_OS_query.h5Seurat", dest = "h5ad", overwrite = T)
## Prep canine PBMC reference
seu.obj <- readRDS("../external_data/GSE225599_final_dataSet_H.rds")
seu.obj$celltype.l3 <- as.character(seu.obj$celltype.l3)
seu.obj <- DietSeurat(seu.obj, assays = "RNA")
SaveH5Seurat(seu.obj, filename = "../external_data/GSE225599_final_dataSet_H.h5Seurat", verbose = TRUE, overwrite = T)
Convert("../external_data/GSE225599_final_dataSet_H.h5Seurat", dest = "h5ad", overwrite = T)
## Prep canine PBMC full dataset -- not used b/c dataset too large to efficiently train
seu.obj <- readRDS("../external_data/GSE225599_final_dataSet_HvO.rds")
seu.obj$celltype.l3 <- as.character(seu.obj$celltype.l3)
seu.obj <- DietSeurat(seu.obj, assays = "RNA")
SaveH5Seurat(seu.obj, filename = "../external_data/GSE225599_final_dataSet_HvO.h5Seurat", verbose = TRUE, overwrite = T)
Convert("../external_data/GSE225599_final_dataSet_HvO.h5Seurat", dest = "h5ad", overwrite = T)
## Prep canine PBMC query dataset
seu.obj <- readRDS("../internal_data/20230505_adj_n5n5_12_5mt_QCfiltered_2500_res0.8_dims45_dist0.3_neigh30_S3.rds")
SaveH5Seurat(seu.obj, filename = "../internal_data/canine_PBMC_query.h5Seurat", verbose = TRUE, overwrite = T)
Convert("../internal_data/canine_PBMC_query.h5Seurat", dest = "h5ad", overwrite = T)
