# seoyeon ahn / scRNA-seq integration - harmony (Seurat v5)
# R_4.3.1

## library
library(Seurat) # Seurat_5.4.0
library(patchwork) # patchwork_1.3.2
library(ggplot2) # ggplot_3.5.2
library(clustree) # clustree_0.5.1

pass="options(future.globals.maxSize = 500 * 1024^3) # 500 Gb RAM
library(Seurat)
library(future)
library(readr)
plan(multisession, workers=72)"

# make directory
dir.create("integration/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("integration/rds", recursive = TRUE, showWarnings = FALSE)

source("utils.R")



### 1. Data Loading
# load cleaned RDS file from previous QC step
rds_files <- list.files("QC/rds/", pattern = "\\.rds$", full.names = TRUE)
obj_list <- list()

for (file in rds_files) {

  sample_name <- gsub(".rds", "", basename(file))

  # read RDS
  tmp_obj <- readRDS(file)

  tmp_obj$sample_id <- sample_name

  # save to list
  obj_list[[sample_name]] <- tmp_obj
}



### 2. Merge all samples into a single Seurat object
merged_obj <- merge(
  x = obj_list[[1]],
  y = obj_list[-1],
  add.cell.ids = names(obj_list)
)

# JoinLayers for initial unintegrated analysis
merged_obj[["RNA"]] <- JoinLayers(merged_obj[["RNA"]])

# check merged object
print("===== Check merged object =====")
merged_obj



### 3. Perform analysis without integration
# standard pipleine to observe batch effect before integration
merged_obj <- NormalizeData(merged_obj)
mereged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj)

## determine dimensionality of the dataset
elbowplots <- visualize_elbow(merged_obj)
  pdf("integration/figures/elbowplot.unintegrated.pdf", width = 10, height = 6)
    print(elbowplots$elbowplot)
    print(elbowplots$elbowplot_adv)
  dev.off()


## clustering without integration for comparison
merged_obj <- FindNeighbors(merged_obj, dims = 1:20, reduction = "pca")
merged_obj <- FindClusters(merged_obj, resolution = 2, cluster.name = "unintegrated_clusters")
merged_obj <- RunUMAP(merged_obj, dims = 1:20, reduction = "pca", reduction.name = "umap.unintegrated")

  dimplot <- DimPlot(merged_obj, reduction = "umap.unintegrated", group.by = c("sample_id", "unintegrated_clusters"))

  pdf("integration/figures/UMAP.unintegrated.pdf", width = 18, height = 6)
    print(dimplot)
  dev.off()



### 4. Perform integration - harmony
## Split layers by sample_id for integration
merged_obj[["RNA"]] <- split(merged_obj[["RNA"]], f = merged_obj$sample_id)

# check splitted object
print("===== Check splitted object =====")
merged_obj

## Run harmony integration
# corrects batch effects while preseving biological variation
merged_obj <- IntegrateLayers(object = merged_obj, 
                              method = HarmonyIntegration, 
                              orig.reduction = "pca", 
                              new.reduction = "harmony", 
                              verbose = FALSE)

# re-join layers after integration
merged_obj[["RNA"]] <- JoinLayers(merged_obj[["RNA"]])



### 5. Integrated Clustering
merged_obj <- FindNeighbors(merged_obj, reduction = "harmony", dims = 1:20)

res_range <- seq(0.1,1,0.1)
merged_obj <- FindClusters(merged_obj, resolution = res_range)
merged_obj <- RunUMAP(merged_obj, dims = 1:20, reduction = "harmony", reduction.name = "umap")

## Save clustering plots by resolution
  pdf("integration/figures/UMAP.integrated.pdf", width = 10, height = 7)
    print(DimPlot(merged_obj, reduction = "umap", group.by = "sample_id"))
  dev.off()

  pdf("integration/figures/UMAP.integrated.clustering.pdf", width = 8, height = 7)
    for(i in res_range) {
      res_col <- paste0("RNA_snn_res.",i)      
      print(DimPlot(merged_obj, reduction = "umap", group.by = res_col, label = TRUE, label.size = 3) + ggtitle(paste0("cluster resolution: ",i)))
    }
  # clustree
    print(clustree(merged_obj, prefix="RNA_snn_res."))
  dev.off()



### 6. Grouping replicates into biological groups and Custom Coloring
sample_ids <- unique(merged_obj$sample_id)

merged_obj$sample_group <- "Unknown"

merged_obj$sample_group[grep("N1", merged_obj$sample_id)] <- "N1"
merged_obj$sample_group[grep("N2", merged_obj$sample_id)] <- "N2"
merged_obj$sample_group[grep("K1_", merged_obj$sample_id)] <- "K1"
merged_obj$sample_group[grep("K1.5", merged_obj$sample_id)] <- "K1.5"
merged_obj$sample_group[grep("K2", merged_obj$sample_id)] <- "K2"
merged_obj$sample_group[grep("K3", merged_obj$sample_id)] <- "K3"
merged_obj$sample_group[grep("K4", merged_obj$sample_id)] <- "K4"
merged_obj$sample_group[grep("K5", merged_obj$sample_id)] <- "K5"
merged_obj$sample_group[grep("K6", merged_obj$sample_id)] <- "K6"

merged_obj$sample_group <- factor(merged_obj$sample_group, levels = c("N1", "N2", "K1", "K1.5", "K2", "K3", "K4", "K5", "K6"))

# assign a custom color palette
# colors are adopted from the reference paper to ensure consistency for comparative analysis
sample_colors <- c(
  "N1" = "#DAF7A6",
  "N2" = "#A1BE6C",
  "K1" = "#EED04D",
  "K1.5" = "#D4AA47",
  "K2" = "#1D89E3",
  "K3" = "#F8BADE",
  "K4" = "#D16BA7",
  "K5" = "#BD14DA",
  "K6" = "#63009B"
)


# save  plots
  pdf("integration/figures/UMAP.final_groups.pdf", width = 8, height = 7)
      print(DimPlot(merged_obj, reduction = "umap", group.by = "sample_group", cols = sample_colors))
  dev.off()


saveRDS(merged_obj, "integration/rds/integration.rds") # save rds

