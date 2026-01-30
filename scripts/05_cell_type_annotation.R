# seoyeon ahn / re-grouping and final cluster annotation
# R_4.3.1

library(Seurat) # Seurat_5.4.0
library(patchwork) # patchwork_1.3.2
library(ggplot2) # ggplot_3.5.2
library(dplyr) # dplyr_1.1.4

pass="options(future.globals.maxSize = 500 * 1024^3) # 500 Gb RAM
library(Seurat)
library(future)
library(readr)
plan(multisession, workers=72)"

# make directory
dir.create("annotation/figures", recursive = TRUE, showWarning=FALSE)
dir.create("annotation/rds", recursive = TRUE, showWarning=FALSE)

source("utils.R")

# read Seurat object
obj <- readRDS("integration/rds/integration.rds")
print(obj)


## Define final annotation marker gene list
marker_list <- list(
  "Acinar" = c("ZG16", "CPA1", "REG1", "PRSS2"),
  "Ductal" = c("KRT19", "SOX9", "CLU", "MMP7", "TSPAN8", "LCN2", "HNF1B"),
  "Nes_Progenitor" = c("NES", "CD44", "VIM", "CDKN2A", "MSN"),
  "Tuft" = c("POU2F3", "DCLK1", "TRPM5"),
  "Neuroendocrine" = c("CHGA", "CHGB", "NEUROD1", "SYP", "PPY", "GCG", "INS1", "INS2", "SST"),
  "Gastric" = c("MUC1", "MUC6", "GKN3", "MUC5AC", "TFF2", "TFF1", "AGR2", "ANXA10"),
  "metaplasia" = c("FOXQ1", "ONECUT2"),
  "PDAC-specific proliferating cell" = c("MKI67", "CDK1"),
  "Fibroblast" = c("ACTA2", "COL1A1", "LUM", "PDGFRB", "COL6A1"),
  "Myeloid cell" = c("NLRP3", "S100A8", "TREM1"),
  "Lymphoid cell" = c("IGHM", "IGLC1", "IGLC3", "MZB1", "IGKC", "CD247", "CD3G", "CD3E", "NKG7", "GIMAP4")
)

all_markers <- unlist(marker_list)



## Re-grouping
# divide the previous cluster 9 into two clusters (early ADM and ADM)
obj$new_manual_groups <- obj$RNA_snn_res.1

Idents(obj) <- "new_manual_groups"

grouping <- list(
  "1" = c(12, 24),
  "2" = c(17, 20),
  "3" = c(14),
  "4" = c(15),
  "5" = c(0, 4, 19),
  "6" = c(1, 3, 7),
  "7" = c(2),
  "8" = c(11),
  "9" = c(13, 18),
  "10" = c(5),
  "11" = c(8, 9, 10, 22),
  "12" = c(6, 16, 23),
  "13" = c(21)
)

new_ids <- setNames(rep(names(grouping), sapply(grouping, length)), unlist(grouping))

obj <- RenameIdents(obj, new_ids)
obj$new_manual_groups <- Idents(obj)




## Visualization
# DotPlot
pdf("annotation/figures/DotPlot.annotation.pdf", width = 18, height = 8)
  print(DotPlot(obj, features = marker_list, dot.scale = 8) + RotatedAxis())
dev.off()

# FeaturePlot
multi_feature_plot(
  object = obj,
  features = all_markers,
  file_name = "annotation/figures/FeaturePlot.annotation.pdf"
  )



## Assigning cell type identity to clusters
new.cluster.ids <- c(
  "Neuroendocrine", 
  "Myeloid cell", 
  "Lymphoid cell", 
  "Tuft cell", 
  "Nes+ progenitor", 
  "Gastric-like cell", 
  "Ductal-like(ADM)", 
  "PDAC-specific proliferating cell", 
  "ADM", 
  "Early ADM", 
  "Acinar", 
  "Fibroblast", 
  "Endothelial cell")

names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)

# final DimPlot
pdf("annotation/figures/DimPlot.final.pdf", width = 10, height = 7)
  print(DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 0.1) + 
        xlab("UMAP 1") + ylab("UMAP 2") +
        theme(axis.title = element_text(size = 18), legend.text = element_text(size = 15)) +
        guides(colour = guide_legend(override.aes = list(size = 3))))
dev.off()

  
saveRDS(obj, "annotation/rds/annotation.rds") # save rds




