# README for Seurat-based Single-Cell RNA Sequencing Analysis

## Description

This R script performs a comprehensive analysis of single-cell RNA sequencing (scRNA-seq) data using the Seurat package. The analysis involves various steps such as loading Seurat objects, visualizing data through dimensionality reduction and plotting, identifying markers for specific cell types, and generating various plots (e.g., feature plots, violin plots) for gene expression.

Key tasks performed by the script:
- Loading and processing Seurat objects.
- Visualization of gene expression with UMAP-based plots.
- Calculation of gene expression percentages across different experimental conditions.
- Identification of differentially expressed markers for cell populations.
- Exporting plots and marker data for further interpretation.

## Prerequisites

Before running the script, make sure to install the necessary R packages. You can do this by running:

```r
install.packages(c("data.table", "readr", "dplyr", "ggplot2"))
install.packages("Seurat")
install.packages("Matrix")
install.packages("openxlsx")
```

## File Structure

- **Input**: 
    - A Seurat object stored in an `.rds` file (`rosmap_seurat.rds`).
    - A predefined list of genes for analysis (e.g., "APOD", "ARFIP2", etc.).
- **Output**:
    - Processed gene data saved as a `.csv` file (`genes_procesados.csv`).
    - Plots for each gene, saved as PNG images in the directories `plots/`, `featureplot_mic/`, and `vlnplot_mic/`.
    - Markers for different cell types (Astrocytes, Endothelial cells, etc.) saved in an Excel file.

## Script Overview

### Step 1: Load and Inspect Seurat Object
The script begins by loading the Seurat object from a `.rds` file and setting the active identity of cells as "seurat_clusters." The object is then inspected using basic functions like `head()` and `str()`.

```r
seurat <- readRDS(file = "rosmap_seurat.rds")
Idents(seurat) <- "seurat_clusters"
head(Idents(seurat))
DimPlot(seurat, reduction = "umap", label = TRUE)
```

### Step 2: Gene Data Processing
A predefined list of genes of interest is created, and their expression data is saved in a CSV file.

```r
genes <- data.frame( ... )  # Predefined list of genes
write_csv(genes, ruta_salida)
```

### Step 3: Gene Expression Plots
For each gene, the script generates:
- Violin plots (`VlnPlot`)
- Feature plots (`FeaturePlot`)

These plots are saved in the `plots/` directory for easy access.

```r
for (gene in genes) {
  vln_plot <- VlnPlot(seurat, features = gene, pt.size = 0.1)
  ggsave(filename = paste0(output_dir, "/", gene, "_vlnplot.png"), plot = vln_plot, width = 6, height = 4)
  
  feature_plot <- FeaturePlot(seurat, features = gene)
  ggsave(filename = paste0(output_dir, "/", gene, "_featureplot.png"), plot = feature_plot, width = 6, height = 4)
}
```

### Step 4: Differential Expression Analysis
The script defines several gene sets corresponding to different experimental conditions (e.g., E2 vs. E3, E4 vs. E2, etc.). It calculates the percentage expression of these genes in each condition and adds it as metadata to the Seurat object.

```r
E2vsE3_percent <- Matrix::colSums(seurat_E2vsE3@assays$RNA@counts) / Matrix::colSums(seurat@assays$RNA@counts)
```

Feature plots for each experimental condition are generated using `FeaturePlot`.

### Step 5: Identifying Markers for Cell Types
The script performs differential expression analysis for several cell types, including microglia, astrocytes, endothelial cells, and more, using the `FindMarkers` function. It then visualizes the results in UMAP plots.

```r
markers_Mic <- FindMarkers(seurat, ident.1 = "Mic", logfc.threshold = 0.25, min.pct = 0.1)
```

For each identified marker, a subset of the data is created, and clustering is performed for further analysis.

### Step 6: Exporting Results
The markers for each cell type are saved into an Excel workbook using the `openxlsx` package. This workbook contains the markers for various cell types, which can be used for downstream analysis.

```r
library(openxlsx)
wb <- createWorkbook()
# Export markers to workbook
```

## Outputs

1. **Gene Expression Plots**: 
    - Feature and violin plots for each gene are saved in the respective directories.
    
2. **Marker Identification**: 
    - A workbook (`markers.xlsx`) containing markers for different cell types (e.g., Microglia, Astrocytes, Endothelial cells, etc.).
    
3. **Cell Type Analysis**:
    - Differentially expressed genes for specific cell types (e.g., Microglia, Astrocytes) are identified using `FindMarkers`.

## File Outputs Example

- `featureplot_mic/FeaturePlot_APOD.png`: Feature plot for the gene "APOD" in microglia cells.
- `vlnplot_mic/VlnPlot_ARFIP2.png`: Violin plot for the gene "ARFIP2" in microglia cells.
- `genes_procesados.csv`: CSV file containing gene expression data.
- `markers.xlsx`: Excel file with markers for different cell types.

## Notes

- The script assumes that the Seurat object is already created and saved as `rosmap_seurat.rds`. This object should contain single-cell RNA-seq data.
- Genes of interest are predefined in the script, but they can be modified based on the user's analysis needs.
- The script generates various types of plots (e.g., UMAP, violin plots) to visualize gene expression patterns across different experimental conditions and cell types.
- Differential gene expression analysis is performed using `FindMarkers`, which compares gene expression between different cell populations.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

---

This README provides a clear outline for running the script, interpreting its outputs, and understanding the analysis workflow.
