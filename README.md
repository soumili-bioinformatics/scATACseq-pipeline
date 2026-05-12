# scATACseq-pipeline

Bioinformatics project for single-cell ATAC-seq data analysis including quality control, dimensionality reduction, clustering, peak analysis, and visualization.

---

## Project Workflow

1. Generate SMAR feature matrix
2. Load Seurat object
3. Fix barcode mismatch
4. Create chromatin assay
5. Perform TF-IDF normalization
6. Run SVD dimensionality reduction
7. Harmony batch correction
8. Clustering analysis
9. Generate UMAP visualization

---

## Tools Used

- R
- Seurat
- Signac
- Harmony
- ggplot2
- dplyr

---

## Repository Structure

```text
scripts/
    run_smar_matrix.R
    harmony_umap_pipeline.R

results/
    umap_sample.png
    umap_clusters.png
    umap_celltype.png
```

---

## Output

- Integrated single-cell ATAC-seq object
- UMAP visualization
- Cluster analysis
- Chromatin accessibility profiling

---

## Author

Soumili Paul

Bioinformatics | Single-cell Omics | Computational Biology
