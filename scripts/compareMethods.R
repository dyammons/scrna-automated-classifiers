#!/usr/bin/env bash

### Preliminaries
library(Seurat)
library(tidyverse)
outName <- "compare"
colz <- list(
  "sc-type" = "#B0CAE1",
  "SingleR" = "#87B13F",
  "Seurat" = "#20A0DA",
  "CellTypist" = "#08B0A1"
)

### Load in results of PBMC classifiers
tissueName <- "PBMC"
## Load in the query dataset
seu.obj <- readRDS("../internal_data/20230505_adj_n5n5_12_5mt_QCfiltered_2500_res0.8_dims45_dist0.3_neigh30_S3.rds")
## Load Seurat predictions
outName <- "Seurat"
annos <- read.csv(paste0("../output/", outName, "/", tissueName, "_anno_predict.csv"))
seu.obj <- AddMetaData(seu.obj, setNames(annos$x, annos$X), paste0(outName, "_predicted.id"))
DimPlot(seu.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = paste0(outName, "_predicted.id"))    
ggsave(paste0("../output/compare/", outName, "_", tissueName, "_anno_UMAP.png"), width = 14, height = 7)
## Load scType predictions
outName <- "scType"
annos <- read.csv(paste0("../output/", outName, "/", tissueName, "_anno_predict.csv"))
seu.obj <- AddMetaData(seu.obj, setNames(annos$x, annos$X), paste0(outName, "_predicted.id"))
DimPlot(seu.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = paste0(outName, "_predicted.id"))    
ggsave(paste0("../output/compare/", outName, "_", tissueName, "_anno_UMAP.png"), width = 14, height = 7)
## Load scType predictions -- highres
annos <- read.csv(paste0("../output/", outName, "/", tissueName, "_highRes_anno_predict.csv"))
seu.obj <- AddMetaData(seu.obj, setNames(annos$x, annos$X), paste0(outName, "_highRes_predicted.id"))
DimPlot(seu.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = paste0(outName, "_highRes_predicted.id"))    
ggsave(paste0("../output/compare/", outName, "_", tissueName, "_highRes_anno_UMAP.png"), width = 14, height = 7)
## Load scType predictions -- highres; refined
annos <- read.csv(paste0("../output/", outName, "/", tissueName, "_highRes_refined_anno_predict.csv"))
seu.obj <- AddMetaData(seu.obj, setNames(annos$x, annos$X), paste0(outName, "_highRes_refined_predicted.id"))
DimPlot(seu.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = paste0(outName, "_highRes_refined_predicted.id"))    
ggsave(paste0("../output/compare/", outName, "_", tissueName, "_highRes_refined_anno_UMAP.png"), width = 14, height = 7)
## Load celltypist predictions
outName <- "celltypist"
annos <- read.csv(paste0("../output/", outName, "/", tissueName, "/", "predicted_labels.csv"))
seu.obj <- AddMetaData(seu.obj, setNames(annos$predicted_labels, annos$X), paste0(outName, "_predicted.id"))
DimPlot(seu.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = paste0(outName, "_predicted.id"))    
ggsave(paste0("../output/compare/", outName, "_", tissueName, "_anno_UMAP.png"), width = 14, height = 7)
## Load singleR predictions
outName <- "singleR"
annos <- read.csv(paste0("../output/", outName, "/", tissueName, "_cluster_anno_predict.csv"))
seu.obj <- AddMetaData(seu.obj, setNames(annos$x, annos$X), paste0(outName, "_cluster_predicted.id"))
DimPlot(seu.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = paste0(outName, "_cluster_predicted.id"))    
ggsave(paste0("../output/compare/", outName, "_", tissueName, "_cluster_anno_UMAP.png"), width = 14, height = 7)
## Load singleR predictions - no Cycling T cells
annos <- read.csv(paste0("../output/", outName, "/", tissueName, "_cluster_noCycle_anno_predict.csv"))
seu.obj <- AddMetaData(seu.obj, setNames(annos$x, annos$X), paste0(outName, "_cluster_noCycle_predicted.id"))
DimPlot(seu.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = paste0(outName, "_cluster_noCycle_predicted.id"))    
ggsave(paste0("../output/compare/", outName, "_", tissueName, "_cluster_noCycle_anno_UMAP.png"), width = 14, height = 7)
## Load singleR predictions - defult cell basis
annos <- read.csv(paste0("../output/", outName, "/", tissueName, "_default_2step_anno_predict.csv"))
seu.obj <- AddMetaData(seu.obj, setNames(annos$x, annos$X), paste0(outName, "_default_2step_predicted.id"))
DimPlot(seu.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = paste0(outName, "_default_2step_predicted.id"))    
ggsave(paste0("../output/compare/", outName, "_", tissueName, "_default_2step_anno_UMAP.png"), width = 14, height = 7)
## Load singleR predictions - assited cell basis
annos <- read.csv(paste0("../output/", outName, "/", tissueName, "_assist_2step_anno_predict.csv"))
seu.obj <- AddMetaData(seu.obj, setNames(annos$x, annos$X), paste0(outName, "_assist_2step_predicted.id"))
DimPlot(seu.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = paste0(outName, "_assist_2step_predicted.id"))    
ggsave(paste0("../output/compare/", outName, "_", tissueName, "_assist_2step_anno_UMAP.png"), width = 14, height = 7)
## Load the "true" cell type annotations that are avalible
outName <- "compare"
seu.anno <- readRDS("../external_data/GSE225599_final_dataSet_HvO.rds")
seu.anno$joinNames <- paste0(substr(colnames(seu.anno), 1, 16), "_", seu.anno$orig.ident)
seu.obj$joinNames <- paste0(substr(colnames(seu.obj), 1, 16), "_", seu.obj$orig.ident)
seu.obj@meta.data <- seu.obj@meta.data %>% 
  rownames_to_column() %>%
  left_join(seu.anno@meta.data[ , c("joinNames", "celltype.l3")], by = "joinNames") %>%
  column_to_rownames()
DimPlot(
  subset(seu.obj, cells = names(seu.obj$celltype.l3[! is.na(seu.obj$celltype.l3)])), 
  reduction = "umap", label = TRUE, repel = TRUE, group.by = "celltype.l3"
) + ggtitle("Ground truth of k9 PBMC query")
ggsave(paste0("../output/", outName, "/", tissueName, "_query_truth_UMAP.png"), width = 14, height = 7)
## Assess classifier accuracy & precision with 5 base approaches
metaData <- seu.obj@meta.data
metaDataZ <- metaData %>% 
  rownames_to_column() %>%
  select(
   celltype.l3, Seurat_predicted.id, scType_predicted.id, celltypist_predicted.id, singleR_cluster_predicted.id
  )
metaDataZ$celltype.l3 <- as.character(metaDataZ$celltype.l3)
metaDataZ <- metaDataZ[! is.na(metaDataZ$celltype.l3), ]
metaDataZ$y <- apply(metaDataZ, 1, function(x) { length(unique(x)) })
table(metaDataZ$y)
metaDataZ$singleR_cluster_predicted.id[is.na(metaDataZ$singleR_cluster_predicted.id)] <- "Unknown" #fix NAs
colnames(metaDataZ)[1] <- "ground_truth"
colz <- c("ground_truth", "Seurat_predicted.id", "scType_predicted.id", "celltypist_predicted.id", "singleR_cluster_predicted.id")
overlap <- unlist(lapply(colz, function(x){
  lapply(colz, function(y){
    sum(metaDataZ[ , x] == metaDataZ[ , y])
  })
}))
mat <- matrix(
  round(overlap / nrow(metaDataZ), 2), 
  ncol = 5, nrow = 5,
  dimnames = list(colz, colz)
)
mat <- mat[rev(rownames(mat)), ]
mat[row(mat) + col(mat) > nrow(mat) + 1] <- NA
melted_cormat <- reshape2::melt(mat, na.rm = TRUE)
#Plot heatmap
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(
   low = "white", high = "red", limit = c(0,1), space = "Lab", 
   name = "% overlapping\nannotations"
 ) +
 theme_minimal() + 
 theme(
   axis.text.x = element_text(angle = 45, vjust = 1, 
   size = 12, hjust = 1)
 ) +
 coord_fixed() + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()
  )
ggsave(paste0("../output/", outName, "/", tissueName, "_coor.png"), width = 7, height = 7)

## Assess classifier accuracy & precision -- with highres
metaData <- seu.obj@meta.data
metaDataZ <- metaData %>% 
  rownames_to_column() %>%
  select(
   celltype.l3, Seurat_predicted.id, scType_predicted.id, scType_highRes_predicted.id, celltypist_predicted.id, singleR_cluster_predicted.id
  )
metaDataZ$celltype.l3 <- as.character(metaDataZ$celltype.l3)
metaDataZ <- metaDataZ[! is.na(metaDataZ$celltype.l3), ]
metaDataZ$y <- apply(metaDataZ, 1, function(x) { length(unique(x)) })
table(metaDataZ$y)
metaDataZ$singleR_cluster_predicted.id[is.na(metaDataZ$singleR_cluster_predicted.id)] <- "Unknown" #fix NAs
colnames(metaDataZ)[1] <- "ground_truth"
colz <- c("ground_truth", "Seurat_predicted.id", "scType_predicted.id", "scType_highRes_predicted.id", "celltypist_predicted.id", "singleR_cluster_predicted.id")
overlap <- unlist(lapply(colz, function(x){
  lapply(colz, function(y){
    sum(metaDataZ[ , x] == metaDataZ[ , y])
  })
}))
mat <- matrix(
  round(overlap / nrow(metaDataZ), 2), 
  ncol = 6, nrow = 6,
  dimnames = list(colz, colz)
)
mat <- mat[rev(rownames(mat)), ]
mat[row(mat) + col(mat) > nrow(mat) + 1] <- NA
melted_cormat <- reshape2::melt(mat, na.rm = TRUE)
# Heatmap
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "white", high = "red",
                      limit = c(0,1), space = "Lab", 
   name = "% overlapping\nannotations") +
  theme_minimal() + 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed() + 
geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank()
)
ggsave(paste0("../output/", outName, "/", tissueName, "_coor_highRes.png"), width = 7, height = 7)

## Assess classifier accuracy & precision -- with highres; refined
metaData <- seu.obj@meta.data
metaDataZ <- metaData %>% 
  rownames_to_column() %>%
  select(
   celltype.l3, Seurat_predicted.id, scType_predicted.id, scType_highRes_predicted.id, scType_highRes_refined_predicted.id, celltypist_predicted.id, singleR_cluster_predicted.id
  )
metaDataZ$celltype.l3 <- as.character(metaDataZ$celltype.l3)
metaDataZ <- metaDataZ[! is.na(metaDataZ$celltype.l3), ]
metaDataZ$y <- apply(metaDataZ, 1, function(x) { length(unique(x)) })
table(metaDataZ$y)
metaDataZ$singleR_cluster_predicted.id[is.na(metaDataZ$singleR_cluster_predicted.id)] <- "Unknown" #fix NAs
colnames(metaDataZ)[1] <- "ground_truth"
colz <- c("ground_truth", "Seurat_predicted.id", "scType_predicted.id", "scType_highRes_predicted.id", "scType_highRes_refined_predicted.id", "celltypist_predicted.id", "singleR_cluster_predicted.id")
overlap <- unlist(lapply(colz, function(x){
  lapply(colz, function(y){
    sum(metaDataZ[ , x] == metaDataZ[ , y])
  })
}))
mat <- matrix(
  round(overlap / nrow(metaDataZ), 2), 
  ncol = 7, nrow = 7,
  dimnames = list(colz, colz)
)
mat <- mat[rev(rownames(mat)), ]
mat[row(mat) + col(mat) > nrow(mat) + 1] <- NA
melted_cormat <- reshape2::melt(mat, na.rm = TRUE)
# Heatmap
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "white", high = "red",
                      limit = c(0,1), space = "Lab", 
   name = "% overlapping\nannotations") +
  theme_minimal() + 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed() + 
geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank()
)
ggsave(paste0("../output/", outName, "/", tissueName, "_coor_highRes_refined.png"), width = 7, height = 7)

res <- lapply(colz[c(2, 3, 6, 7)], function(x){
  metaDataZ <- metaDataZ[ , c(colz[1], x)]
  colnames(metaDataZ)[2] <- "predicted"
  pred.df <- lapply(unique(metaDataZ$ground_truth), function(y){
    res <- metaDataZ %>% 
      summarize(
        TP = sum(ground_truth == y & predicted == y),
        TN = sum(ground_truth != y & predicted != y),
        FP = sum(ground_truth != y & predicted == y),
        FN = sum(ground_truth == y & predicted != y),
        predictor = x,
        celltype = y
      )
    rownames(res) <- paste0(x, "_", y)
    return(res)
  })
  res <- do.call(rbind, pred.df)
  res %>% 
    mutate(
      recall = TP / (TP + FP),
      precision = TP / (TP + FN),
      F1 = round(2 * (precision * recall) / (precision + recall), 2)
    )
})
res <- do.call(rbind, res)
## Make the plot
ggplot(data = res, aes(predictor, celltype, fill = F1))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "white", high = "red",
                      limit = c(0,1), space = "Lab", 
   name = "F1") +
  theme_minimal() + 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed() + 
geom_text(aes(predictor, celltype, label = F1), color = "black", size = 4) +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank()
) + coord_fixed(ratio = 0.2)
ggsave(paste0("../output/", outName, "/", tissueName, "_f1.png"), width = 10, height = 10)

lapply(colz[2:5], function(predictor){
  classyHelper.df <- as.data.frame(table(metaDataZ[ , "ground_truth"], metaDataZ[ , predictor])) %>%
    group_by(Var1) %>%
    mutate(
      totalFreq = sum(Freq),
      Percent = round(Freq / totalFreq, 2)
    ) %>%
    ungroup()
  classyHelper.df$Percent[classyHelper.df$Freq == 0] <- NA
  ## Make the plot
  ggplot(data = classyHelper.df, aes(Var1, Var2, fill = Percent))+
   geom_tile(color = "white") +
    viridis::scale_fill_viridis() + 
  #  scale_fill_gradient2(low = "white", high = "red",
  #                       limit = c(0,1), space = "Lab", 
  #    name = "F1") +
    theme_minimal() + 
    labs(x = "Ground truth", y = "Predicted cell type", title = predictor) + 
   theme(axis.text.x = element_text(angle = 45, vjust = 1, 
      size = 12, hjust = 1))+
   coord_fixed() + 
  # geom_text(aes(Var1, Var2, label = Freq), color = "black", size = 4) +
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()
  )
  ggsave(paste0("../output/", outName, "/", tissueName, "_", predictor, "_predHelp.png"), width = 10, height = 10)
})

## Replot with more linent assessment of accuracy
metaData <- seu.obj@meta.data
metaDataZ <- metaData %>% 
  rownames_to_column() %>%
  select(
   celltype.l3, Seurat_predicted.id, scType_predicted.id, celltypist_predicted.id, singleR_predicted.id
  )
metaDataZ$celltype.l3 <- as.character(metaDataZ$celltype.l3)
metaDataZ <- metaDataZ[! is.na(metaDataZ$celltype.l3), ]
metaDataZ$y <- apply(metaDataZ, 1, function(x) { length(unique(x)) })
table(metaDataZ$y)
metaDataZ$singleR_predicted.id[is.na(metaDataZ$singleR_predicted.id)] <- "Unknown" #fix NAs
colnames(metaDataZ)[1] <- "ground_truth"
metaDataZ <- as.data.frame(sapply(metaDataZ, function(x){
  x[grepl("TEM|TCM", x)] <- "T_mem"  
  return(x)
}))
colz <- c("ground_truth", "Seurat_predicted.id", "scType_predicted.id", "celltypist_predicted.id", "singleR_predicted.id")
overlap <- unlist(lapply(colz, function(x){
  lapply(colz, function(y){
    sum(metaDataZ[ , x] == metaDataZ[ , y])
  })
}))
mat <- matrix(
  round(overlap / nrow(metaDataZ), 2), 
  ncol = 5, nrow = 5,
  dimnames = list(colz, colz)
)
mat <- mat[rev(rownames(mat)), ]
mat[row(mat) + col(mat) > nrow(mat) + 1] <- NA
melted_cormat <- reshape2::melt(mat, na.rm = TRUE)
# Heatmap
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "white", high = "red",
                      limit = c(0,1), space = "Lab", 
   name = "% overlapping\nannotations") +
  theme_minimal() + 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed() + 
geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank()
)
ggsave(paste0("../output/", outName, "/", tissueName, "_coor_revisited.png"), width = 7, height = 7)


### Load in results of cross-tissue classifiers
## Load in the query dataset
seu.obj <- readRDS("../internal_data/20230504_naiveTX_n6_n6_QCfiltered_1_2500_res0.8_dims45_dist0.3_neigh30_S3.rds")
## Load celltypist predictions and generate final plot
tissueName <- "OS_PBMC"
outName <- "celltypist"
annos <- read.csv(paste0("../output/", outName, "/", tissueName, "/", "predicted_labels.csv"))
seu.obj <- AddMetaData(seu.obj, setNames(annos$predicted_labels, annos$X), paste0(outName, "_predicted.id"))
DimPlot(seu.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = paste0(outName, "_predicted.id"))    
ggsave(paste0("../output/", outName, "/", tissueName, "/", tissueName, "_anno_UMAP.png"), width = 14, height = 7)

tissueName <- "OS_PBMC_prob_match"
outName <- "celltypist"
annos <- read.csv(paste0("../output/", outName, "/", tissueName, "/", "predicted_labels.csv"))
seu.obj <- AddMetaData(seu.obj, setNames(annos$predicted_labels, annos$X), paste0(outName, "_predicted.id"))
seu.obj$celltypist_predicted.id <- as.factor(
  ifelse(
    grepl("Unassigned|\\|", seu.obj$celltypist_predicted.id), 
    NA, 
    seu.obj$celltypist_predicted.id
  )
)
DimPlot(seu.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = paste0(outName, "_predicted.id"))    
ggsave(paste0("../output/", outName, "/", tissueName, "/", tissueName, "_anno_UMAP.png"), width = 14, height = 7)
## Load in the T cell subcluster query dataset
seu.obj <- readRDS("../internal_data/harmony_Tcells_res0.6_dims30_dist0.3_neigh30_S3.rds")
## Load seurat predictions and generate final plot
tissueName <- "OStcell_PBMC"
outName <- "Seurat"
annos <- read.csv(paste0("../output/", outName, "/", tissueName, "_anno_predict.csv"))
seu.obj <- AddMetaData(seu.obj, setNames(annos$x, annos$X), paste0(outName, "_predicted.id"))
DimPlot(seu.obj, reduction = "umap.integrated.harmony", label = TRUE, repel = TRUE, group.by = paste0(outName, "_predicted.id"))
ggsave(paste0("../output/", outName, "/", tissueName, "_anno_UMAP.png"), width = 10, height = 7)
DimPlot(seu.obj, reduction = "umap.integrated.harmony", label = TRUE, repel = TRUE, group.by = "clusterID2_integrated.harmony") + NoLegend()
ggsave(paste0("../output/", outName, "/", tissueName, "_cluster_UMAP.png"), width = 7, height = 7)
seu.obj$ground_truth <- seu.obj$celltype.l3_naive
DimPlot(seu.obj, reduction = "umap.integrated.harmony", label = TRUE, repel = TRUE, group.by = "ground_truth")
ggsave(paste0("../output/", outName, "/", tissueName, "_cluster_UMAP.png"), width = 10, height = 7)

classyHelper.df <- as.data.frame(table(seu.obj$ground_truth, seu.obj$Seurat_predicted.id)) %>%
  group_by(Var1) %>%
  mutate(
    totalFreq = sum(Freq),
    Percent = round(Freq / totalFreq, 2)
  ) %>%
  ungroup()
classyHelper.df$Percent[classyHelper.df$Freq == 0] <- NA
## Make the plot
ggplot(data = classyHelper.df, aes(Var1, Var2, fill = Percent))+
 geom_tile(color = "white") +
  viridis::scale_fill_viridis() + 
  theme_minimal() + 
  labs(x = "Ground truth", y = "Predicted cell type", title = "Seurat_predicted.id") + 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed() + 
theme(
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank()
)
ggsave(paste0("../output/", outName, "/", tissueName, "_Seurat_predicted.id_predHelp.png"), width = 10, height = 10)



## Load Seurat predictions
outName <- "Seurat"
annos <- read.csv(paste0("../output/", outName, "/", tissueName, "_anno_predict.csv"))
seu.obj <- AddMetaData(seu.obj, setNames(annos$x, annos$X), paste0(outName, "_predicted.id"))
## Load scType predictions
outName <- "scType"
annos <- read.csv(paste0("../output/", outName, "/", tissueName, "_anno_predict.csv"))
seu.obj <- AddMetaData(seu.obj, setNames(annos$x, annos$X), paste0(outName, "_predicted.id"))
## Load celltypist predictions and generate final plot
outName <- "celltypist"
annos <- read.csv(paste0("../output/", outName, "/", tissueName, "/", "predicted_labels.csv"))
seu.obj <- AddMetaData(seu.obj, setNames(annos$predicted_labels, annos$X), paste0(outName, "_predicted.id"))
DimPlot(seu.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = paste0(outName, "_predicted.id"))    
ggsave(paste0("../output/", outName, "/", tissueName, "/", tissueName, "_anno_UMAP.png"), width = 14, height = 7)

outName <- "compare"
metaData <- seu.obj@meta.data

scType_names <- c(
  "conventional DC subtype 2", "CD8 exhausted", "CD4 na ve", 
  "Endothelial cell", "Malignant osteoblast subtype 3", "Malignant osteoblast subtype 2", 
  "Cycling T cell", "Hypoxic osteoblast", "Lipid associated TAM (SPP2 high)", 
  "Neutrophil", "Malignant osteoblast subtype 1", "Lipid associated TAM (C1QC high)", 
  "Fibroblast", "T-IFN", "Cycling osteoblast subtype 1", "Intermediate TAM", 
  "Mast cell", "CD4- tumor infiltrating monocyte", "Mature osteoclast", 
  "IFN-osteoblast", "plasmacytoid DC", "mature regulatory DC"
)
corrected_names <- c(
  "cDC2", "CD8_ex", "CD4_naive", "Endothelial cell", "Osteoblast_3", 
  "Osteoblast_2", "T_cycling", "Hypoxic_osteoblast", "LA-TAM_SPP2_hi", 
  "Neutrophil", "Osteoblast_1", "LA-TAM_C1QC_hi", "Fibroblast", "T_IFN", 
  "Osteoblast_cycling_1", "TAM_INT", "Mast cell", "CD4-_TIM", "Mature_OC", 
  "IFN-osteoblast", "pDC", "mregDC"
)

namez <- setNames(scType_names, corrected_names)
matches <- lapply(namez, function(x){which(metaData$scType_predicted.id == x)})
for (i in 1:(length(matches))) {
  metaData$scType_predicted.id[unlist(unname(matches[i]))] <- names(matches[i])
}

metaDataZ <- metaData %>% 
  select(
   celltype.l3_naive2, Seurat_predicted.id, scType_predicted.id, celltypist_predicted.id
  ) %>%
  filter(celltype.l3_naive2 != "NA")

metaDataZ$y <- apply(metaDataZ, 1, function(x) { length(unique(x)) })
table(metaDataZ$y)

colz <- c("celltype.l3_naive2", "Seurat_predicted.id", "scType_predicted.id", "celltypist_predicted.id")
overlap <- unlist(lapply(colz, function(x){
  lapply(colz, function(y){
    sum(metaDataZ[ , x] == metaDataZ[ , y])
  })
}))

mat <- matrix(
  round(overlap / nrow(metaDataZ), 2), 
  ncol = 4, nrow = 4,
  dimnames = list(colz, colz)
)

melted_cormat <- reshape2::melt(mat, na.rm = TRUE)
# Heatmap
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "white", high = "red",
                      limit = c(0,1), space = "Lab", 
   name = "Accuracy (%)") +
  theme_minimal() + 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed() + 
geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank()
)
ggsave(paste0("../output/", outName, "/", tissueName, "_coor.png"), width = 7, height = 7)



