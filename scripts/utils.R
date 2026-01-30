# seoyeon ahn / single cell utils
# R_4.3.1

## library
library(Seurat) # Seurat_5.4.0
library(patchwork) # patchwork_1.3.2
library(ggplot2) # ggplot2_3.5.2


## Load single cell data and creat Seurat object
load_sc_data <- function(sample_path, sample_name) {

  barcodes_csv <- list.files(sample_path, pattern = "barcodes", full.names = TRUE)
  genes_csv <- list.files(sample_path, pattern = "genes", full.names = TRUE)
  mtx_file <- list.files(sample_path, pattern = "molecule_counts.mtx.gz", full.names = TRUE)

  barcodes <- read.csv(gzfile(barcodes_csv), header = FALSE, stringsAsFactors = FALSE, colClasses = c("NULL", "character"))[, 1]
  genes <- read.csv(gzfile(genes_csv), header = FALSE, stringsAsFactors = FALSE, colClasses = c("NULL", "character"))[, 1]
  genes <- make.unique(genes)

  mtx <- readMM(gzfile(mtx_file))
  mtx <- as(mtx, "dgCMatrix") # cells x genes

  counts <- t(mtx)
  rownames(counts) <- genes
  colnames(counts) <- barcodes

  # record raw information
  qc_df1[sample_name,"features"] <<- dim(counts)[1]
  qc_df1[sample_name,"cells"] <<- dim(counts)[2]

  seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 3, min.features = 200)

  return(seurat_obj)
}


## Doublet Finding and Removing using scDblFinder
doublet_clean <- function(seurat_obj, sample_name) {

  # convert to SCE for scDblFinder
  sce <- SingleCellExperiment(list(counts = GetAssayData(seurat_obj, layer = "counts")))
  sce <- scDblFinder(sce, dbr = 0.075)

  seurat_obj$scDblFinder_class <- sce$scDblFinder.class
  seurat_obj$scDblFinder_score <- sce$scDblFinder.score
  
  # standard pipeline for doublet visualization (temporary object)
  seu <- NormalizeData(seurat_obj)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  seu <- ScaleData(seu)

  seu <- RunPCA(seu)
  seu <- RunUMAP(seu, dims = 1:10)

  # doublet plot
  n_singlet <- sum(seu$scDblFinder_class == "singlet")
  n_doublet <- sum(seu$scDblFinder_class == "doublet")


  plot_labels <- c(paste0("Singlet (n=", n_singlet, ")"),
                   paste0("Doublet (n=", n_doublet, ")"))
  names(plot_labels) <- c("singlet", "doublet")

  dimplot <- DimPlot(seu, group.by = "scDblFinder_class",
                     cols = c("singlet"="grey", "doublet"="red")) +
             scale_color_manual(values = c("singlet" = "grey", "doublet" = "red"),
                                labels = plot_labels) +
            ggtitle(paste0(sample_name, " Doublet Detection"))

  pdf(paste0("QC/qc/DoubletPlot.", sample_name, ".pdf"), width = 7, height = 7)
    print(dimplot)
  dev.off()

  # excluding expected doublets
  seurat_obj_clean <- subset(seu, subset = scDblFinder_class == "singlet")

  return(seurat_obj_clean)
}



## Custom Elbow Plot
visualize_elbow <- function(object) {
  plot1 <- ElbowPlot(object)

  # elbow detection logic
  pct <- object[["pca"]]@stdev / sum(object[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  pcs <- min(co1, co2)

  plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))
  plot2 <- ggplot(plot_df, aes(cumu, pct, label = rank, color = rank  > pcs)) +
    geom_text() +
    geom_vline(xintercept = 90, color = "grey") +
    geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
    theme_bw()

  return(list(elbowplot = plot1, elbowplot_adv = plot2))
}


## Custom multi-page FeaturePlot for large gene lists
multi_feature_plot <- function(object, features, file_name, n_per_page = 9, ncol = 3) {

  # calculate number of pages needed
  feature_split <- split(features, ceiling(seq_along(features) / n_per_page))
  nrow_val <- n_per_page / ncol

  pdf(file_name, width = 12, height = 12)

    for (i in seq_along(feature_split)) {
      current_features <- feature_split[[i]]
      n_features <- length(current_features)

      # generate plot for the current subset of genes
      plot_list <- FeaturePlot(
        object = object,
        features = current_features,
        reduction = "umap",
        pt.size = 0.2,
        order = T,
        combine = FALSE
      )

      if (n_features < n_per_page) {
        n_spacers <- n_per_page - n_features
        for (s in 1:n_spacers) {
          plot_list[[n_features + s]] <- plot_spacer()
        }
      }

      combined_plot <- wrap_plots(plot_list, ncol = ncol, nrow = nrow_val) +
                       plot_annotation(theme = theme(plot.title = element_text(size = 15)))
      print(combined_plot)
    }
  dev.off()
}

