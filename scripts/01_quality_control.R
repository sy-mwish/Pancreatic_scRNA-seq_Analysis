# seoyeon ahn / quality control, doublet removal (scDblFinder)
# R_4.3.1

## library
library(Seurat) # Seurat_5.4.0
library(Matrix) # Matrix_1.6.5
library(ggplot2) # ggplot2_3.5.2
library(scDblFinder) #scDblFinder_1.16.0
library(clustree) # clustree_0.5.1
library(dplyr) # dplyr_1.1.4
library(cowplot) # cowplot_1.2.0
library(SingleCellExperiment) # SingleCellExperiment_1.24.0

pass="options(future.globals.maxSize = 500 * 1024^3) # 500 Gb RAM
library(Seurat)
library(future)
library(readr)
plan(multisession, workers=72)"


sample_dir <- "/DATA2/seoyeon/Epithelial"
# make directory
dir.create("QC/qc/Scaled_QC",recursive = TRUE, showWarning=FALSE)
dir.create("QC/rds", recursive = TRUE, showWarning=FALSE)

source("utils.R")


### Main analysis loop
smpls <- list.files(sample_dir)

# QC Table
qc_df1 <- data.frame(features = numeric(length(smpls)), cells = numeric(length(smpls)), row.names = smpls)
qc_df2 <- data.frame(features = numeric(length(smpls)), cells = numeric(length(smpls)), row.names = smpls)
qc_df3 <- data.frame(features = numeric(length(smpls)), cells = numeric(length(smpls)), row.names = smpls)


for (smpl in smpls) {
  print(paste("Basic Analysis for", smpl, "is processing..."))
  
  sample_path <- file.path(sample_dir, smpl)
  
  ## 1. Load and initial QC
  seurat_obj <- load_sc_data(sample_path, smpl)

  seurat_obj$percent.mt <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

  # rawQC plotting
  vlnplot <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") + theme(legend.position="none")
  plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  fsplot <- plot1 + plot2

  pdf(paste0("QC/qc/RawQC.feature_count_percentMT.",smpl,".pdf"), width = 10, height = 7)
    print(vlnplot)
    print(fsplot)
  dev.off()


  ## 2. Filtering
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 500 &
                                            nCount_RNA > 1500 &
                                            percent.mt < 15)
  
  # record information after filtering
  qc_df2[smpl,"features"] <- dim(seurat_obj)[1]
  qc_df2[smpl,"cells"] <- dim(seurat_obj)[2]


  ## 3. Doublet Cleaning
  seurat_obj <- doublet_clean(seurat_obj, smpl)

  # record final information
  qc_df3[smpl,"features"] <- dim(seurat_obj)[1]
  qc_df3[smpl,"cells"] <- dim(seurat_obj)[2]
  
  
  ## 4. Scaling & PCA
  seurat_obj <- ScaleData(seurat_obj, rownames(seurat_obj), verbose = FALSE)
  # perform linear dimension reduction
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), verbose = FALSE)
    vizdimplot <- VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca")
    dimplot <- DimPlot(seurat_obj, reduction = "pca") + NoLegend()
    elbowplot <- ElbowPlot(seurat_obj)
  
  pdf(paste("QC/qc/Scaled_QC/dimensional_reduction." , smpl, ".pdf"), width = 8, height = 7)
    print(vizdimplot)
    print(dimplot)
    print(elbowplot)
  dev.off()

  pdf(paste("QC/qc/Scaled_QC/dimensional_reduction.DimHeatmap.",smpl,".pdf"), width = 18, height = 15)
    DimHeatmap(seurat_obj, dims = 1, cells = 500, balanced = TRUE)
    DimHeatmap(seurat_obj, dims = 1:9, ncol=3, cells = 500, balanced = TRUE)
  dev.off()

  ## 5. Clustering the cells
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:10)
  seurat_obj <- FindClusters(seurat_obj, resolution = seq(0.1,1,0.1))
  # run non-linear dimensional reduction (UMAP/tSNE)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)


  # save clustering plots by resolution
  pdf(paste0("QC/qc/UMAP.RNA_snn_resolution.",smpl,".pdf"), width = 10, height = 8)
    for(i in seq(0.1,1,0.1)){
      print(DimPlot(seurat_obj, reduction = "umap", group.by = paste0("RNA_snn_res.",i)))
    }
  # clustree
    print(clustree(seurat_obj, prefix="RNA_snn_res."))
  dev.off()


  ## 6. Save Results
  saveRDS(seurat_obj, paste0("QC/rds/", smpl, ".rds"))
}

## QC Table
qc_table <- cbind(qc_df1, qc_df2, qc_df3)
print(qc_table)

# save 
write.csv(qc_table, "QC/qc/qc_table.csv", row.names = TRUE)


