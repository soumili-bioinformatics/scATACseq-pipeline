library(Signac)
library(Seurat)
library(GenomicRanges)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(harmony)
library(ggplot2)
library(dplyr)

brca <- readRDS("BRCA.rds")
smar.matrix <- readRDS("smar_matrix.rds")

# Barcode correction
seurat.cells <- colnames(brca)
seurat.barcodes <- sub(".*_", "", seurat.cells)

common.barcodes <- intersect(colnames(smar.matrix), seurat.barcodes)

smar.matrix <- smar.matrix[, common.barcodes]

match.index <- match(common.barcodes, seurat.barcodes)
colnames(smar.matrix) <- seurat.cells[match.index]

# Add SMAR assay
annotation <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"

brca[["SMAR"]] <- CreateChromatinAssay(
  counts = smar.matrix,
  fragments = Fragments(brca),
  annotation = annotation
)

# Normalization
DefaultAssay(brca) <- "SMAR"

brca <- RunTFIDF(brca)
brca <- FindTopFeatures(brca, min.cutoff = "q0")
brca <- RunSVD(brca)

# Harmony integration
harmony_embeddings <- harmony::HarmonyMatrix(
  data_mat = Embeddings(brca, "lsi")[,2:20],
  meta_data = brca@meta.data,
  vars_use = "Piece_ID",
  do_pca = FALSE
)

brca[["harmony"]] <- CreateDimReducObject(
  embeddings = harmony_embeddings,
  key = "harmony_",
  assay = DefaultAssay(brca)
)

# Clustering
brca <- FindNeighbors(
  brca,
  reduction = "harmony",
  dims = 1:19
)

brca <- FindClusters(brca, resolution = 0.5)

# UMAP
brca <- RunUMAP(
  brca,
  reduction = "harmony",
  dims = 1:19
)

# Plotting
p1 <- DimPlot(
  brca,
  reduction = "umap",
  group.by = "Piece_ID",
  pt.size = 0.1
)

p2 <- DimPlot(
  brca,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  pt.size = 0.1
)

p3 <- DimPlot(
  brca,
  reduction = "umap",
  group.by = "cell_type",
  shuffle = TRUE,
  label = TRUE,
  pt.size = 0.1
)

# Save plots
ggsave("umap_sample.png", p1)
ggsave("umap_clusters.png", p2)
ggsave("umap_celltype.png", p3)

# Save object
saveRDS(brca, "brca_after_harmony_umap.rds")
