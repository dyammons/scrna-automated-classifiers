#!/usr/bin/env Rscript

#Load packages and sc-type functions
library(Seurat)
library(HGNChelper)
library(tidyverse)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
outName <- "scType"

### PBMC analysis
## Run sc-type cluster based approach
#Load in the data
seu.obj <- readRDS("../internal_data/20230505_adj_n5n5_12_5mt_QCfiltered_2500_res0.8_dims45_dist0.3_neigh30_S3.rds")
scRNAseqData_scaled <- as.matrix(seu.obj@assays$integrated$scale.data)
tissueName <- "PBMC"
#Format cell type gene signatures - only positive ct markers
gene_sig <- read.csv("../external_data/supplemental_data_4.csv")
modulez <- split(gene_sig$gene, gene_sig$cellType_l2)
gs_list <- list(
  "gs_positive" = modulez,
  "gs_negative" =  lapply(modulez, function(mod){mod <- character(0)})
)
#Run sctype on the integrated data matrix
es.max <- sctype_score(
  scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative
)
#Organize the results by cluster
cL_resutls <- do.call("rbind", lapply(unique(seu.obj@meta.data$clusterID), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seu.obj@meta.data[seu.obj@meta.data$clusterID==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seu.obj@meta.data$clusterID==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
#Set low ScType scored clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
#Update Seurat object metadata
Idents(seu.obj) <- "clusterID"
seu.obj <- RenameIdents(seu.obj, setNames(sctype_scores$type, sctype_scores$cluster))
seu.obj$sctype_classification <- Idents(seu.obj)
#Plot the results
DimPlot(seu.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification')    
ggsave(paste0("../output/", outName, "/", tissueName, "_anno_UMAP.png"), width = 14, height = 7)
write.csv(seu.obj$sctype_classification, file = paste0("../output/", outName, "/", tissueName, "_anno_predict.csv"))
## Run sc-type cluster based approach increase resolution for demonstration purposes
tissueName <- "PBMC_highRes"
#Load data and increase clustering resolution
seu.obj <- readRDS("../internal_data/20230505_adj_n5n5_12_5mt_QCfiltered_2500_res0.8_dims45_dist0.3_neigh30_S3.rds")
seu.obj <- FindClusters(seu.obj, resolution = 1.6, algorithm = 3, graph.name = "integrated_snn", cluster.name = "clusterID_highRes")
#Plot high res clustering results
DimPlot(seu.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = "clusterID_highRes") + 
  ggtitle("Unsupervised clustering (high res) of k9 PBMC query") + 
  NoLegend()
ggsave(paste0("../output/", outName, "/", tissueName, "_query_UMAP_highRes.png"), width = 7, height = 7)
scRNAseqData_scaled <- as.matrix(seu.obj@assays$integrated$scale.data)
#Format cell type gene signatures - only positive ct markers
gene_sig <- read.csv("../external_data/supplemental_data_4.csv")
modulez <- split(gene_sig$gene, gene_sig$cellType_l2)
gs_list <- list(
    "gs_positive" = modulez,
    "gs_negative" =  lapply(modulez, function(mod){mod <- character(0)})
    )
#Run sctype on the integrated data matrix
es.max <- sctype_score(
  scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative
)
# Organize the results by cluster
cL_resutls <- do.call("rbind", lapply(unique(seu.obj@meta.data$clusterID_highRes), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seu.obj@meta.data[seu.obj@meta.data$clusterID_highRes==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seu.obj@meta.data$clusterID_highRes==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
#Set low ScType scored clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
#Update Seurat object metadata
Idents(seu.obj) <- "clusterID_highRes"
seu.obj <- RenameIdents(seu.obj, setNames(sctype_scores$type, sctype_scores$cluster))
seu.obj$sctype_classification <- Idents(seu.obj)
#Plot the results
DimPlot(seu.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification')    
ggsave(paste0("../output/", outName, "/", tissueName, "_anno_UMAP.png"), width = 14, height = 7)
write.csv(seu.obj$sctype_classification, file = paste0("../output/", outName, "/", tissueName, "_anno_predict.csv"))
## Run sc-type cluster based approach increase resolution for demonstration purposes
tissueName <- "PBMC_highRes_refined"
#Load data and increase clustering resolution as done above
seu.obj <- readRDS("../internal_data/20230505_adj_n5n5_12_5mt_QCfiltered_2500_res0.8_dims45_dist0.3_neigh30_S3.rds")
seu.obj <- FindClusters(seu.obj, resolution = 1.6, algorithm = 3, graph.name = "integrated_snn", cluster.name = "clusterID_highRes")
scRNAseqData_scaled <- as.matrix(seu.obj@assays$integrated$scale.data)
#Format cell type gene signatures - add in positive and negative ct markers
gene_sig <- read.csv("../external_data/supplemental_data_4.csv")
modulez <- split(gene_sig$gene, gene_sig$cellType_l2)
gs_list <- list(
    "gs_positive" = modulez,
    "gs_negative" =  lapply(modulez, function(mod){mod <- character(0)})
    )
gs_list$gs_positive$`CD4+ Naive` <- c(gs_list$gs_positive$`CD4+ Naive`, "CD4")
gs_list$gs_positive$`CD8+ Naive` <- c(gs_list$gs_positive$`CD8+ Naive`, "CD8A")
gs_list$gs_negative$`CD4+ Naive` <- c("CD8A")
gs_list$gs_negative$`CD8+ Naive` <- c("CD4")
gs_list$gs_negative$`NK cell` <- c("CD3E")
#Run sctype on the integrated data matrix
es.max <- sctype_score(
  scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative
)
#Organize the results by cluster
cL_resutls <- do.call("rbind", lapply(unique(seu.obj@meta.data$clusterID_highRes), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seu.obj@meta.data[seu.obj@meta.data$clusterID_highRes==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seu.obj@meta.data$clusterID_highRes==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
#Set low ScType scored clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
#Update Seurat object metadata
Idents(seu.obj) <- "clusterID_highRes"
seu.obj <- RenameIdents(seu.obj, setNames(sctype_scores$type, sctype_scores$cluster))
seu.obj$sctype_classification <- Idents(seu.obj)
#Plot the results
DimPlot(seu.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification')    
ggsave(paste0("../output/", outName, "/", tissueName, "_anno_UMAP.png"), width = 14, height = 7)
write.csv(seu.obj$sctype_classification, file = paste0("../output/", outName, "/", tissueName, "_anno_predict.csv"))

### PBMC on OS cross tissue analysis
tissueName <- "OS_PBMC"
#Load in the data
seu.obj <- readRDS("../internal_data/20230504_naiveTX_n6_n6_QCfiltered_1_2500_res0.8_dims45_dist0.3_neigh30_S3.rds")
scRNAseqData_scaled <- as.matrix(seu.obj@assays$integrated$scale.data)
#Format ct gene signatures
gene_sig <- read.csv("../external_data/supplemental_data_4.csv")
modulez <- split(gene_sig$gene, gene_sig$cellType_l2)
gs_list <- list(
  "gs_positive" = modulez,
  "gs_negative" =  lapply(modulez, function(mod){mod <- character(0)})
)
#Run sctype on the integrated data matrix
es.max <- sctype_score(
  scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative
)
#Organize the results by cluster
cL_resutls <- do.call("rbind", lapply(unique(seu.obj@meta.data$clusterID), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seu.obj@meta.data[seu.obj@meta.data$clusterID==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seu.obj@meta.data$clusterID==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
#Set low ScType scored clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells] <- "Unknown"
#Update Seurat object metadata
Idents(seu.obj) <- "clusterID"
seu.obj <- RenameIdents(seu.obj, setNames(sctype_scores$type, sctype_scores$cluster))
seu.obj$sctype_classification <- Idents(seu.obj)
#Plot the results
DimPlot(seu.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification')    
ggsave(paste0("../output/", outName, "/", tissueName, "_anno_UMAP.png"), width = 14, height = 7)
write.csv(seu.obj$sctype_classification, file = paste0("../output/", outName, "/", tissueName, "_anno_predict.csv"))

### BONUS OS analysis
tissueName <- "OS"
#Load in the data
seu.obj <- readRDS("../internal_data/20230504_naiveTX_n6_n6_QCfiltered_1_2500_res0.8_dims45_dist0.3_neigh30_S3.rds")
scRNAseqData_scaled <- as.matrix(seu.obj@assays$integrated$scale.data)
#Format ct gene signatures
gene_sig <- read.csv("../external_data/supplemental_data_2_OS_atlas.csv")
modulez <- split(gene_sig$gene, gene_sig$celltype.l3)
gs_list <- list(
    "gs_positive" = modulez,
    "gs_negative" =  lapply(modulez, function(mod){mod <- character(0)})
    )
#Run sctype on the integrated data matrix
es.max <- sctype_score(
  scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative
)
#Organize the results by cluster
cL_resutls <- do.call("rbind", lapply(unique(seu.obj@meta.data$clusterID), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seu.obj@meta.data[seu.obj@meta.data$clusterID==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seu.obj@meta.data$clusterID==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
#Set low ScType scored clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
#Update Seurat object metadata
Idents(seu.obj) <- "clusterID"
seu.obj <- RenameIdents(seu.obj, setNames(sctype_scores$type, sctype_scores$cluster))
seu.obj$sctype_classification <- Idents(seu.obj)
#Plot the results
DimPlot(seu.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification')    
ggsave(paste0("../output/", outName, "/", tissueName, "_anno_UMAP.png"), width = 14, height = 7)
write.csv(seu.obj$sctype_classification, file = paste0("../output/", outName, "/", tissueName, "_anno_predict.csv"))
