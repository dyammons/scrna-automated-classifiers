#!/usr/bin/env bash

# Load packages and sc-type functions
library(Seurat)
library(SingleR)
library(tidyverse)
outName <- "singleR"

### PBMC analysis
#Prep the data
reference <- readRDS("../external_data/GSE225599_final_dataSet_H.rds")
seu.obj <- readRDS("../internal_data/20230505_adj_n5n5_12_5mt_QCfiltered_2500_res0.8_dims45_dist0.3_neigh30_S3.rds")
cntData <- GetAssayData(seu.obj, layer = "data", assay = "RNA")
refData <- GetAssayData(reference, layer = "data", assay = "RNA")
cntData <- cntData[intersect(rownames(cntData), rownames(refData)), ]
refData <- refData[intersect(rownames(cntData), rownames(refData)), ]
## Run SingleR cluster based approach
tissueName <- "PBMC_cluster"
rst <- SingleR(test = cntData, ref = refData, clusters = seu.obj$clusterID, labels = reference$celltype.l3)
seu.obj$singleR_classification <- rst$pruned.labels[match(seu.obj@meta.data[["clusterID"]], rownames(rst))]
p <- DimPlot(seu.obj, reduction = "umap", group.by = "singleR_classification", label = TRUE)
ggsave(paste0("../output/", outName, "/", tissueName, "_anno_UMAP.png"), width = 14, height = 7)
write.csv(seu.obj$singleR_classification, file = paste0("../output/", outName, "/", tissueName, "_anno_predict.csv"))
## Run SingleR cluster based approach, but with "Cycling T cell" removed
tissueName <- "PBMC_cluster_noCycle"
reference <- subset(reference, invert = T, cells = names(reference$celltype.l3[reference$celltype.l3 == "Cycling T cell"]))
reference$celltype.l3 <- droplevels(reference$celltype.l3)
refData <- GetAssayData(reference, layer = "data", assay = "RNA")
refData <- refData[intersect(rownames(cntData), rownames(refData)), ]
rst <- SingleR(test = cntData, ref = refData, clusters = seu.obj$clusterID, labels = reference$celltype.l3)
seu.obj$singleR_classification <- rst$pruned.labels[match(seu.obj@meta.data[["clusterID"]], rownames(rst))]
p <- DimPlot(seu.obj, reduction = "umap", group.by = "singleR_classification", label = TRUE)
ggsave(paste0("../output/", outName, "/", tissueName, "_anno_UMAP.png"), width = 14, height = 7)
write.csv(seu.obj$singleR_classification, file = paste0("../output/", outName, "/", tissueName, "_anno_predict.csv"))
## Two step approach with default settings
#Prep the data
tissueName <- "PBMC_default_2step"
reference <- readRDS("../external_data/GSE225599_final_dataSet_H.rds")
seu.obj <- readRDS("../internal_data/20230505_adj_n5n5_12_5mt_QCfiltered_2500_res0.8_dims45_dist0.3_neigh30_S3.rds")
cntData <- GetAssayData(seu.obj, layer = "data", assay = "RNA")
refData <- GetAssayData(reference, layer = "data", assay = "RNA")
cntData <- cntData[intersect(rownames(cntData), rownames(refData)), ]
refData <- refData[intersect(rownames(cntData), rownames(refData)), ]
#Train
ref_trained <- trainSingleR(
  refData,
  reference$celltype.l3
)
#Classify
rst <- classifySingleR(
  test = cntData,
  trained = ref_trained
)
#Plot
seu.obj$singleR_classification <- rst$pruned.labels[match(colnames(seu.obj), rownames(rst))]
p <- DimPlot(seu.obj, reduction = "umap", group.by = "singleR_classification", label = TRUE)
ggsave(paste0("../output/", outName, "/", tissueName, "_anno_UMAP.png"), width = 14, height = 7)
write.csv(seu.obj$singleR_classification, file = paste0("../output/", outName, "/", tissueName, "_anno_predict.csv"))
## Two step approach with cell type markers
#Prep the data
tissueName <- "PBMC_assist_2step"
reference <- readRDS("../external_data/GSE225599_final_dataSet_H.rds")
seu.obj <- readRDS("../internal_data/20230505_adj_n5n5_12_5mt_QCfiltered_2500_res0.8_dims45_dist0.3_neigh30_S3.rds")
cntData <- GetAssayData(seu.obj, layer = "data", assay = "RNA")
refData <- GetAssayData(reference, layer = "data", assay = "RNA")
cntData <- cntData[intersect(rownames(cntData), rownames(refData)), ]
refData <- refData[intersect(rownames(cntData), rownames(refData)), ]
#Load cell type gene signatures 
gene_sig <- read.csv("../external_data/supplemental_data_4.csv")
gene_sig <- split(gene_sig$gene, gene_sig$cellType_l2)
#Train
ref_trained <- trainSingleR(
  refData,
  reference$celltype.l3,
  genes = gene_sig
)
#Classify
rst <- classifySingleR(
  test = cntData,
  trained = ref_trained
)
#Plot
seu.obj$singleR_classification <- rst$pruned.labels[match(colnames(seu.obj), rownames(rst))]
p <- DimPlot(seu.obj, reduction = "umap", group.by = "singleR_classification", label = TRUE)
ggsave(paste0("../output/", outName, "/", tissueName, "_anno_UMAP.png"), width = 14, height = 7)
write.csv(seu.obj$singleR_classification, file = paste0("../output/", outName, "/", tissueName, "_anno_predict.csv"))

### BONUS OS analysis
#Prep the data
seu.obj <- readRDS("../internal_data/20230504_naiveTX_n6_n6_QCfiltered_1_2500_res0.8_dims45_dist0.3_neigh30_S3.rds")
reference <- readRDS("../external_data/canine_naive_n6_annotated.rds")
tissueName <- "OS_cluster"
cntData <- GetAssayData(seu.obj, layer = "data", assay = "RNA")
refData <- GetAssayData(reference, layer = "data", assay = "RNA")
cntData <- cntData[intersect(rownames(cntData), rownames(refData)), ]
refData <- refData[intersect(rownames(cntData), rownames(refData)), ]
## Run SingleR cluster based approach
rst <- SingleR(test = cntData, ref = refData, clusters = seu.obj$clusterID, labels = reference$celltype.l3)
seu.obj$singleR_classification <- rst$pruned.labels[match(seu.obj@meta.data[["clusterID"]], rownames(rst))]
p <- DimPlot(seu.obj, reduction = "umap", group.by = "singleR_classification", label = TRUE)
ggsave(paste0("../output/", outName, "/", tissueName, "_anno_UMAP.png"), width = 14, height = 7)
write.csv(seu.obj$singleR_classification, file = paste0("../output/", outName, "/", tissueName, "_anno_predict.csv"))
## Two step approach with cell type markers
#Prep the data
seu.obj <- readRDS("../internal_data/20230504_naiveTX_n6_n6_QCfiltered_1_2500_res0.8_dims45_dist0.3_neigh30_S3.rds")
reference <- readRDS("../external_data/canine_naive_n6_annotated.rds")
tissueName <- "OS_assist_2step"
cntData <- GetAssayData(seu.obj, layer = "data", assay = "RNA")
refData <- GetAssayData(reference, layer = "data", assay = "RNA")
cntData <- cntData[intersect(rownames(cntData), rownames(refData)), ]
refData <- refData[intersect(rownames(cntData), rownames(refData)), ]
#Load cell type gene signatures and correct naming discprency
gene_sig <- read.csv("../external_data/supplemental_data_2_OS_atlas.csv")
orig.name <- c(
  "ANGIO-TAM", "Activated TAM", "B cell", "CD320 osteoclast", 
  "CD4 activated", "CD4 follicular helper", "CD4 na ve", "CD4 regulatory", 
  "CD4+ tumor infiltrating monocyte", "CD4- tumor infiltrating monocyte", 
  "CD8 SPP1+", "CD8 effector", "CD8 exhausted", "Cycling T cell", 
  "Cycling osteoblast subtype 1", "Cycling osteoblast subtype 2", 
  "Cycling osteoblast subtype 3", "Cycling osteoblast subtype 4", 
  "Cycling osteoclast 1/2", "Endothelial cell", "Fibroblast", "Hypoxic osteoblast", 
  "IFN-TAM", "IFN-osteoblast", "Intermediate TAM", "Lipid associated TAM (C1QC high)", 
  "Lipid associated TAM (SPP2 high)", "Malignant osteoblast subtype 1", 
  "Malignant osteoblast subtype 2", "Malignant osteoblast subtype 3", 
  "Mast cell", "Mature osteoclast", "Neutrophil", "Plasma cell", 
  "T-IFN", "conventional DC subtype 1", "conventional DC subtype 2", 
  "mature regulatory DC", "plasmacytoid DC", "precursor DC"
 )
new.name <- c(
  "ANGIO_TAM", "TAM_ACT", "B cell", "CD320_OC", "CD4_act", "CD4_fh", "CD4_naive",
  "CD4_reg", "CD4+_TIM", "CD4-_TIM", "CD8_SPP1_hi", "CD8_eff", "CD8_ex",
  "T_cycling", "Osteoblast_cycling_1", "Osteoblast_cycling_2", "Osteoblast_cycling_3", 
  "Osteoblast_cycling_4", "Cycling_OC", "Endothelial cell", "Fibroblast", 
  "Hypoxic_osteoblast", "IFN-TAM", "IFN-osteoblast",  "TAM_INT", "LA-TAM_C1QC_hi", 
  "LA-TAM_SPP2_hi", "Osteoblast_1", "Osteoblast_2", "Osteoblast_3", "Mast cell", 
  "Mature_OC", "Neutrophil", "Plasma cell", "T_IFN", "cDC1", "cDC2", "mregDC", 
  "pDC", "preDC"
)
namez <- setNames(orig.name, new.name)
matches <- lapply(namez, function(x){which(gene_sig$celltype.l3 == x)})
for (i in 1:(length(matches))) {
  gene_sig$celltype.l3[unlist(unname(matches[i]))] <- names(matches[i])
}
gene_sig <- c(
    split(gene_sig$gene, gene_sig$celltype.l3), 
    list("NK" = c("KLRF1", "STMN2", "PAX4", "NCR3", "F2RL3", "CD96", "IL2RB", "IGSF3", "FREM1", "FASLG"))
)
#Train
ref_trained <- trainSingleR(
  refData,
  reference$celltype.l3,
  genes = gene_sig
)
#Classify
rst <- classifySingleR(
  test = cntData,
  trained = ref_trained
)
#Plot
seu.obj$singleR_classification <- rst$pruned.labels[match(colnames(seu.obj), rownames(rst))]
p <- DimPlot(seu.obj, reduction = "umap", group.by = "singleR_classification", label = TRUE)
ggsave(paste0("../output/", outName, "/", tissueName, "_anno_UMAP.png"), width = 14, height = 7)
write.csv(seu.obj$singleR_classification, file = paste0("../output/", outName, "/", tissueName, "_anno_predict.csv"))
