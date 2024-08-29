#!/usr/bin/env Rscript

# Load libraries
library(Seurat)
library(tidyverse)
outName <- "general"


### Plot PBMC
tissueName <- "PBMC"
## Load in and plot the PBMC query dataset
seu.obj <- readRDS("../internal_data/20230505_adj_n5n5_12_5mt_QCfiltered_2500_res0.8_dims45_dist0.3_neigh30_S3.rds")
DimPlot(seu.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = "clusterID") + 
  ggtitle("Unsupervised clustering of k9 PBMC query") + 
  NoLegend()
ggsave(paste0("../output/", outName, "/", tissueName, "_query_UMAP.png"), width = 7, height = 7)
DimPlot(seu.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = "clusterID") + 
  ggtitle("Unsupervised clustering of k9 PBMC query") + 
  NoLegend()
ggsave(paste0("../output/", outName, "/", tissueName, "_query_UMAP.png"), width = 7, height = 7)
## Load in and plot the PBMC annotated reference dataset
seu.obj <- readRDS("../external_data/GSE225599_final_dataSet_H.rds")
DimPlot(seu.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = "celltype.l3") + 
  ggtitle("Annotated k9 PBMC atlas") + 
  NoLegend()
ggsave(paste0("../output/", outName, "/", tissueName, "_ref_UMAP.png"), width = 7, height = 7)


### Plot OS
tissueName <- "OS"
## Load in and plot the OS query dataset
seu.obj <- readRDS("../internal_data/20230504_naiveTX_n6_n6_QCfiltered_1_2500_res0.8_dims45_dist0.3_neigh30_S3.rds")
DimPlot(seu.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = "clusterID") + 
  ggtitle("Unsupervised clustering of k9 OS query") + 
  NoLegend()
ggsave(paste0("../output/", outName, "/", tissueName, "_query_UMAP.png"), width = 7, height = 7)
FeaturePlot(seu.obj, features = c("PTPRC", "CD3E", "CD4", "CD8A"))
ggsave(paste0("../output/", outName, "/", tissueName, "_query_feats.png"), width = 8, height = 7)
## Load in and plot the OS annotated reference dataset
seu.obj <- readRDS("../external_data/canine_naive_n6_annotated.rds")
DimPlot(seu.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = "celltype.l3") + 
  ggtitle("Annotated k9 OS atlas") + 
  NoLegend()
ggsave(paste0("../output/", outName, "/", tissueName, "_ref_UMAP.png"), width = 7, height = 7)
