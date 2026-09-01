suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(yaml)
})

params <- yaml::read_yaml("configs/analysis_params.yaml")
samples <- readr::read_csv("configs/samples.csv", show_col_types = FALSE)

read_sample <- function(sample_id, group, batch, data_path, notes) {
  counts <- Read10X(data.dir = data_path)
  object <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = params$quality_control$min_cells_per_gene,
    min.features = params$quality_control$min_genes_per_cell
  )
  object$sample_id <- sample_id
  object$group <- group
  object$batch <- batch
  object
}

objects <- purrr::pmap(samples, read_sample)
sc <- merge(objects[[1]], y = objects[-1], add.cell.ids = samples$sample_id)

sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")

qc_plot <- VlnPlot(
  sc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample_id",
  ncol = 3,
  pt.size = 0
)
ggsave("results/figures/seurat_qc_violin.png", qc_plot, width = 12, height = 5)

sc <- subset(
  sc,
  subset =
    nFeature_RNA >= params$quality_control$min_genes_per_cell &
    nFeature_RNA <= params$quality_control$max_genes_per_cell &
    percent.mt <= params$quality_control$max_mito_percent
)

sc <- NormalizeData(sc, scale.factor = params$normalization$scale_factor)
sc <- FindVariableFeatures(sc, nfeatures = params$normalization$n_highly_variable_genes)
sc <- ScaleData(sc)
sc <- RunPCA(sc)
sc <- FindNeighbors(sc, dims = 1:params$dimension_reduction$n_pcs)
sc <- FindClusters(sc, resolution = params$dimension_reduction$clustering_resolution)
sc <- RunUMAP(sc, dims = 1:params$dimension_reduction$n_pcs)

umap_plot <- DimPlot(sc, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  DimPlot(sc, reduction = "umap", group.by = "sample_id")
ggsave("results/figures/seurat_umap_clusters_samples.png", umap_plot, width = 12, height = 5)

markers <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
readr::write_csv(markers, "results/tables/seurat_cluster_markers.csv")

saveRDS(sc, "results/objects/seurat_processed.rds")
