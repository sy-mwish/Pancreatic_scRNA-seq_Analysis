# seoyeon ahn / Grouping clusters based on marker gene expression (start from resolution = 1.0)
# R_4.3.1

library(Seurat) # Seurat_5.4.0
library(patchwork) # patchwork_1.3.2
library(ggplot2) # ggplot_3.5.2

pass="options(future.globals.maxSize = 500 * 1024^3) # 500 Gb RAM
library(Seurat)
library(future)
library(readr)
plan(multisession, workers=72)"

# make directory
dir.create("grouping/figures", recursive = TRUE, showWarning=FALSE)
dir.create("grouping/rds", recursive = TRUE, showWarning=FALSE)

source("utils.R")

# read Seurat object
obj <- readRDS("integration/rds/integration.rds")
print(obj)



## Define marker gene list
marker_list <- list(
  "Acinar" = c("ZG16", "CPA1"),
  "Ductal" = c("KRT19", "SOX9", "CLU", "MMP7", "TSPAN8", "LCN2"),
  "Nes_Progenitor" = c("NES", "CD44", "VIM", "CDKN2A", "MSN"),
  "Tuft" = c("POU2F3", "DCLK1", "TRPM5"),
  "Neuroendocrine" = c("CHGA", "CHGB", "NEUROD1", "SYP", "PPY", "GCG", "INS1", "INS2", "SST"),
  "Gastric" = c("MUC1", "MUC6", "GKN3", "MUC5AC", "TFF2", "TFF1", "AGR2", "ANXA10"),
  "fibrolast" = c("ACTA2", "COL1A1", "LUM"),
  "ADM" = c("FOXQ1", "ONECUT2")
)

all_markers <- unlist(marker_list, use.names = FALSE)



## Check marker genes expression
# FeaturePlot
multi_feature_plot(
  object = obj, 
  features = all_markers, 
  file_name = "grouping/figures/FeaturePlot.markers.pdf", 
  n_per_page = 9, 
  ncol = 3
)



## Cluster Grouping 
# initializing groups from resolution 1.0
obj$manual_groups <- obj$RNA_snn_res.1

Idents(obj) <- "manual_groups"


# mapping table for merging clusters based on expression similarity
grouping <- list(
  "1" = c(12, 24),
  "2" = c(17, 20),
  "3" = c(14),
  "4" = c(15),
  "5" = c(0, 4, 19),
  "6" = c(1, 3, 7),
  "7" = c(2),
  "8" = c(11),
  "9" = c(5, 13, 18),
  "10" = c(8, 9, 10, 22),
  "11" = c(6, 16, 23),
  "12" = c(21)
)


# convert list to named vector for RenameIdents
new_ids <- setNames(rep(names(grouping), sapply(grouping, length)), unlist(grouping))

# re-assigning cluster identities
obj <- RenameIdents(obj, new_ids)
obj$manual_groups <- Idents(obj)



## Visualization
# UMAP of grouped clusters
pdf("grouping/figures/UMAP.grouping.pdf", width = 8, height = 7)
  print(DimPlot(obj, reduction = "umap", group.by = "manual_groups", label = TRUE) + ggtitle("Grouped Clusters"))
dev.off()

# DotPlot & ViolinPlot
dotplot <- DotPlot(obj, features = marker_list, dot.scale = 8) + RotatedAxis()
violinplot <- VlnPlot(obj, features = all_markers, pt.size = 0, combine = FALSE)
wrap_plots(plots = violinplot, ncol = 4)

pdf("grouping/figures/grouping_marker.pdf", width = 14, height = 8)
  print(dotplot)
  print(violinplot)
dev.off()

saveRDS(obj, "grouping/rds/grouping.rds") # save rds




