library(Signac)
library(Seurat)
library(rtracklayer)
library(Matrix)

cat("Loading SMAR regions...\n")
smars <- rtracklayer::import("hg38-all-smars.bed.gz")

cat("Finding fragment files...\n")
frag.files <- list.files(
  "filtered_BRCA",
  pattern = "sorted.tsv.gz$",
  full.names = TRUE
)

cat("Creating fragment objects...\n")
frag.objects <- lapply(
  frag.files,
  function(x) CreateFragmentObject(path = x)
)

cat("Running FeatureMatrix...\n")
smar.matrix <- FeatureMatrix(
  fragments = frag.objects,
  features = smars
)

cat("Saving matrix...\n")
saveRDS(smar.matrix, "smar_matrix.rds")

cat("Finished successfully\n")
