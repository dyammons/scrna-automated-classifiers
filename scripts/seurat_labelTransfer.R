#!/usr/bin/env 

# Load libraries
source("https://raw.githubusercontent.com/dyammons/scrna-seq/main/analysis-code/customFunctions_Seuratv5.R") #custom plotting functions
library(Seurat)
library(tidyverse)
outName <- "Seurat"

## Define function
runSeuratLabTransfer <- function(
  seu.obj,
  reference,
  reduction = "umap",
  outName,
  tissueName
){
  ref.anchors <- FindTransferAnchors(
    reference = reference, query = seu.obj, dims = 1:30, reference.reduction = "pca", 
    features = rownames(seu.obj)[rownames(seu.obj) %in% rownames(reference)]
  )
  predictions <- TransferData(
    anchorset = ref.anchors, 
    refdata = reference$celltype.l3, 
    dims = 1:30
  )
  seu.obj <- AddMetaData(seu.obj, metadata = predictions)
  pi <- DimPlot(
    seu.obj, reduction = reduction, group.by = "predicted.id",
    pt.size = 0.1, label = T, label.box = T, repel = T
  )
  ggsave(paste0("../output/", outName, "/", tissueName, "_anno_UMAP.png"), width = 14, height = 7)
  write.csv(seu.obj$predicted.id, file = paste0("../output/", outName, "/", tissueName, "_anno_predict.csv"))
}


## Run PBMC
seu.obj <- readRDS("../internal_data/20230505_adj_n5n5_12_5mt_QCfiltered_2500_res0.8_dims45_dist0.3_neigh30_S3.rds")
reference <- readRDS("../external_data/GSE225599_final_dataSet_H.rds")
runSeuratLabTransfer(
  seu.obj = seu.obj, reference = reference, reduction = "umap",
  outName = outName, tissueName = "PBMC"
)
## Run OS to annotate with PBMC
seu.obj <- readRDS("../internal_data/20230504_naiveTX_n6_n6_QCfiltered_1_2500_res0.8_dims45_dist0.3_neigh30_S3.rds")
reference <- readRDS("../external_data/GSE225599_final_dataSet_H.rds")
runSeuratLabTransfer(
  seu.obj = seu.obj, reference = reference, reduction = "umap",
  outName = outName, tissueName = "OS_PBMC"
)
## Run OS T cells to annotate with PBMC
seu.obj <- readRDS("../internal_data/20230504_naiveTX_n6_n6_QCfiltered_1_2500_res0.8_dims45_dist0.3_neigh30_S3.rds")
#integrate the data
seu.obj <- subset(seu.obj, subset = clusterID %in% c(7, 9, 26, 28))
seu.obj <- integrateData(
  seu.obj = seu.obj, dout = "../internal_data/", 
  outName = outName, normalization.method = "LogNormalize",
  runAllMethods = FALSE, method = "HarmonyIntegration", 
  new.reduction.name = "integrated.harmony", indReClus = T, 
  vars.to.regress = "percent.mt", saveRDS = F
)
#complete data visualization
seu.obj <- dataVisUMAP(
  seu.obj = seu.obj, outDir = "../internal_data/", outName = "harmony_Tcells", 
  final.dims = 30, final.res = 0.6, stashID = "clusterID2", algorithm = 3, min.dist = 0.3, n.neighbors = 30,
  prefix = "RNA_snn_res.", assay = "RNA", reduction = "integrated.harmony",
  saveRDS = T, return_obj = T, returnFeats = T,
  features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                "IL7R", "FLT3", "S100A12", "CD68", 
                "CD4", "MS4A1", "TOP2A","IL23R")
)
reference <- readRDS("../external_data/GSE225599_final_dataSet_H.rds")
tissueName <- "OStcell_PBMC"
runSeuratLabTransfer(
  seu.obj = seu.obj, reference = reference, reduction = "umap.integrated.harmony",
  outName = outName, tissueName = "OStcell_PBMC"
)
## Load data - tumor test
seu.obj <- readRDS("../internal_data/20230504_naiveTX_n6_n6_QCfiltered_1_2500_res0.8_dims45_dist0.3_neigh30_S3.rds")
reference <- readRDS("../external_data/canine_naive_n6_annotated.rds")
runSeuratLabTransfer(
  seu.obj = seu.obj, reference = reference, reduction = "umap",
  outName = outName, tissueName = "OS"
)