# seoyeon ahn / trajectory analysis 
# root = normal(N1) sample cell
# R_4.3.1

## library
library(monocle3) # monocle3_1.4.26
library(ggplot2) # ggplot2_3.5.2
library(dplyr) # dplyr_1.1.4
library(Seurat) # Seurat_5.4.0
library(SeuratWrappers) # SeuratWrappers_0.4.0

# make directory
dir.create("trajectory/n1_root/figures", recursive = TRUE, showWarnings=FALSE)
dir.create("trajectory/n1_root/cds", recursive = TRUE, showWarnings=FALSE)


## Convert seurat object to monocle CDS object
# read seurat object
obj <- readRDS("grouping/rds/grouping.rds")
print(obj)

cds <- as.cell_data_set(obj, assay = "RNA")

# injecting seurat's UMAP coordinate into CDS
cds@int_colData$reducedDims$UMAP <- obj@reductions$umap@cell.embeddings



## Clustering : assign cluster and partition to each cells
cds <- cluster_cells(cds, reduction_method = "UMAP")
plot <- plot_cells(cds, color_cells_by = "partition")

# Replace  cluster information with seurat's
cds@clusters$UMAP$clusters <- obj@active.ident



## Learn the trajectory graph
# use seurat's coordinate as they are
cds <- learn_graph(cds, use_partition = TRUE)
plot1 <- plot_cells(cds,
                      color_cells_by = "cluster",
                      label_groups_by_cluster=FALSE,
                      label_leaves=FALSE,
                      label_branch_points=FALSE)

plot2 <- plot_cells(cds,
                  color_cells_by = "sample_group",
                  label_cell_groups=FALSE,
                  label_leaves=TRUE,
                  label_branch_points=TRUE,
                  graph_label_size=1.5)

pdf("trajectory/n1_root/figures/learn_graph.pdf", width = 8, height = 6)
  print(plot1)
  print(plot2)
dev.off()



## Order the cells in pseudotime
# select root node
# root node : node with the most cell in the normal(N1) sample
get_root_nodes <- function(cds, sample_col = "sample_group", root_sample = "N1") {

  all_partitions <- unique(partitions(cds))
  root_nodes <- c()

  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])


  for (p in all_partitions) {
    cells_in_p <- names(partitions(cds)[partitions(cds) == p])

    n1_cells_in_p <- cells_in_p[colData(cds)[cells_in_p, sample_col] == root_sample]

    if (length(n1_cells_in_p) > 0) {
      node <- igraph::V(principal_graph(cds)[["UMAP"]])$name[
        as.numeric(names(which.max(table(closest_vertex[n1_cells_in_p, ]))))
      ]
      root_nodes <- c(root_nodes, node)
    }
  }

  root_nodes
}

cds <- order_cells(cds, root_pr_nodes = get_root_nodes(cds))


plot1 <- plot_cells(cds,
                  color_cells_by = "pseudotime",
                  label_cell_groups=FALSE,
                  label_leaves=FALSE,
                  label_branch_points=FALSE,
                  graph_label_size=1.5)

plot2 <- plot_cells(cds,
                  color_cells_by = "pseudotime",
                  label_cell_groups=FALSE,
                  label_leaves=FALSE,
                  label_branch_points=FALSE,
                  graph_label_size=1.5,
                  show_trajectory_graph=FALSE)


pdf("trajectory/n1_root/figures/order_cells.pdf", width = 8, height = 6)
  print(plot1)
  print(plot2)
dev.off()


saveRDS(obj, "trajectory/n1_root/cds/n1_root.rds") # save cds




