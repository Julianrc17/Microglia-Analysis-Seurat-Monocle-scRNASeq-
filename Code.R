# Loading Libraries

library(data.table)
library(Matrix)
library(Seurat)
library(readr)
library(Seurat)
library(dplyr)
library(ggplot2) 

# Loading .rds seurat file (first, create a seurat object)

seurat <-readRDS(file = "rosmap_seurat.rds")
Idents(seurat)<-"seurat_clusters"
head(Idents(seurat))

DimPlot(seurat, reduction = "umap", label = TRUE)

str(seurat)


ruta_salida <- "/home/julian/scRNAseq/genes_procesados.csv"

genes <- data.frame(
  Symbol = c("APOD", "ARFIP2", "ARL1", "ARL2", "ATE1", "ATG7", "BCDIN3D", "BIRC2", "BTG1", "C14orf93",
             "C5orf38", "CABP7", "CCDC24", "CCL25", "CDK5RAP3", "CDKN2B", "CEP20", "CHCHD7", "CLEC4G", "CLUL1",
             "COASY", "COMTD1", "COPS8", "CRP", "CRYZL1", "CTF1", "DCUN1D5", "DHRS9", "DMKN", "EFCAB1",
             "EMILIN3", "FAM20A", "FBLIM1", "FCHSD1", "FOXO1", "GCLM", "GGT2", "HAX1", "HBQ1", "HES1",
             "HSD17B8", "HYKK", "IRF6", "ITPK1", "KCNRG", "KCTD2", "KRT19", "LCN1", "LMO4", "LRRN1",
             "MAGEA3", "MAP1LC3B2", "MAPK6", "MAPK8", "MOB4", "MTMR7", "NFIA", "OTULIN", "PICK1", "PLEKHO2",
             "PSMB10", "RDH16", "RND1", "RNF122", "S100A13", "SAR1A", "SELENOS", "SLC27A2", "SNRPF", "SPATC1L",
             "SPC25", "ST8SIA1", "STAR", "STK10", "TBCA", "TIMM13", "TMCC3", "TP53I11", "TRIM54", "UNG",
             "VPS29", "WNT10B", "ZNF483"),
  
  E2vsE3 = c(NA, NA, NA, -0.81, NA, NA, -0.44, 0.19, -0.27, NA,
             -0.69, NA, 1.17, NA, NA, NA, 0.38, NA, 0.18, NA,
             -0.22, NA, NA, NA, NA, NA, NA, 0.59, 0.27, NA,
             NA, 0.20, NA, NA, 0.33, NA, 0.36, NA, NA, 0.20,
             NA, NA, NA, NA, 0.25, NA, NA, 0.20, NA, NA,
             -0.31, NA, NA, NA, NA, 0.21, NA, NA, NA, NA,
             -0.21, 0.21, NA, NA, -0.93, 0.18, NA, 0.37, NA, 0.60,
             NA, 0.24, -0.27, 0.63, NA, NA, NA, NA, 0.25, 0.94,
             NA, NA, NA),
  
  E4vsE2 = c(NA, -0.19, -0.27, NA, -0.47, NA, NA, -0.22, 0.50, -0.31,
             1.38, NA, -1.72, 0.34, -0.32, -0.20, -0.45, 0.52, -0.20, -0.31,
             NA, -0.23, -0.24, -0.59, -0.27, 1.08, -0.50, -0.76, -0.40, -0.16,
             NA, -0.20, -0.22, -0.25, -0.98, -0.25, -0.25, 0.63, -0.61, -0.21,
             NA, -0.33, -0.21, -0.26, -0.25, NA, -0.21, NA, -0.53, 1.78,
             NA, -0.25, -0.38, 0.18, NA, -0.59, 0.29, 0.28, -0.27, 0.37,
             NA, -0.52, 0.24, 0.21, -0.75, -0.19, -0.60, -0.39, -0.30, -0.64,
             2.48, -0.50, 0.51, -0.74, NA, -0.26, -0.25, -1.27, -0.33, -0.94,
             NA, NA, NA),
  
  E4vsE3 = c(-0.20, NA, NA, -0.91, -0.36, -0.17, -0.51, NA, 0.24, -0.25,
             0.75, -0.20, -0.55, 0.39, NA, NA, NA, 0.41, NA, -0.35,
             -0.21, -0.19, NA, -0.40, -0.25, 0.99, -0.34, NA, NA, NA,
             -0.17, NA, NA, NA, -0.62, NA, NA, 0.54, -0.74, NA,
             -0.19, -0.20, -0.29, -0.29, NA, -0.39, -0.26, NA, -0.54, 1.79,
             -0.36, -0.22, -0.24, NA, -0.31, -0.37, NA, 0.26, NA, 0.24,
             NA, -0.30, 0.25, NA, -1.68, NA, -0.72, NA, -0.42, NA,
             2.44, -0.25, 0.26, NA, -0.99, NA, -0.24, -1.13, NA, NA,
             -0.34, -0.29, -0.19)
)

write_csv(genes, ruta_salida)

genes <- read_csv(ruta_salida)

library(ggplot2)
genes <- c("APOD", "ARFIP2", "ARL1", "ARL2", "ATE1", "ATG7", "BCDIN3D", "BIRC2", "BTG1", "C14orf93",
           "C5orf38", "CABP7", "CCDC24", "CCL25", "CDK5RAP3", "CDKN2B", "CHCHD7", "CLEC4G", "CLUL1",
           "COASY", "COMTD1", "COPS8", "CTF1", "DCUN1D5", "DHRS9", "DMKN", "EFCAB1",
           "EMILIN3", "FAM20A", "FBLIM1", "FCHSD1", "FOXO1", "GCLM", "GGT2", "HAX1", "HBQ1", "HES1",
           "HSD17B8", "HYKK", "IRF6", "ITPK1", "KCNRG", "KCTD2", "KRT19", "LCN1", "LMO4", "LRRN1",
           "MAP1LC3B2", "MAPK6", "MAPK8", "MOB4", "MTMR7", "NFIA", "OTULIN", "PICK1",
           "PSMB10", "RDH16", "RND1", "RNF122", "S100A13", "SAR1A", "SLC27A2", "SNRPF", "SPATC1L",
           "SPC25", "ST8SIA1", "STAR", "STK10", "TBCA", "TIMM13", "TMCC3", "TP53I11", "TRIM54", "UNG",
           "VPS29", "WNT10B", "ZNF483")


output_dir <- "plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

for (gene in genes) {
  vln_plot <- VlnPlot(seurat, features = gene, pt.size = 0.1) + ggtitle(paste("VlnPlot -", gene))
  ggsave(filename = paste0(output_dir, "/", gene, "_vlnplot.png"), plot = vln_plot, width = 6, height = 4)
  
  feature_plot <- FeaturePlot(seurat, features = gene) + ggtitle(paste("FeaturePlot -", gene))
  ggsave(filename = paste0(output_dir, "/", gene, "_featureplot.png"), plot = feature_plot, width = 6, height = 4)
}


genes_problemas <- c("CEP20", "CRP", "CRYZL1", "MAGEA3", "PLEKHO2", "SELENOS")

# Candidates Plot

E2vsE3 <- c("ARL2", "BCDIN3D", "BIRC2", "C5orf38", "CCDC24", "CEP20", "CLEC4G","COASY", "DHRS9", "DMKN", "FAM20A", "FOXO1", "GGT2", "HES1", "KCNRG", "LCN1", "MAGEA3", "MTMR7", "RDH16", "S100A13", "SAR1A", "SLC27A2", "SPATC1L", "ST8SIA1", "STAR", "STK10", "TRIM54", "UNG")
E4vsE2 <- c("ARFIP2", "ARL1", "ATE1", "BIRC2", "BTG1", "C14orf93", "C5orf38", "CCDC24","CCL25","CDK5RAP3","CDKN2B","CEP20","CHCHD7","CLEC4G","CLUL1","COMTD1","COPS8","CRP","CRYZL1","CTF1","DCUN1D5","DHRS9", "DMKN","EFCAB1","FAM20A","FBLIM1","FCHSD1","FOXO1","GCLM","GGT2", "HAX1","HBQ1","HES1","HYKK","IRF6","ITPK1","KCNRG","KRT19","LMO4","LRRN1","MAP1LC3B2","MAPK6","MAPK8","MTMR7","NFIA","OTULIN","PICK1","PLEKHO2","RDH16","RND1","RNF122","S100A13","SAR1A","SELENOS","SLC27A2","SNRPF","SPATC1L","SPC25","ST8SIA1","STAR","STK10","TIMM13","TMCC3","TP53I11","TRIM54","UNG")
E4vsE3 <- c("APOD","ARL2","ATE1","ATG7","BCDIN3D","BTG1","C14orf93","C5orf38","CABP7","CCDC24","CCL25","CHCHD7","CLUL1","COASY","COMTD1","CRP","CRYZL1","CTF1","DCUN1D5","EMILIN3","FOXO1","HAX1","HBQ1","HSD17B8","HYKK","IRF6","ITPK1","KCTD2","KRT19","LMO4","LRRN1","MAGEA3","MAP1LC3B2","MAPK6","MOB4","MTMR7","OTULIN","PLEKHO2","RDH16","RND1","S100A13","SELENOS","SNRPF","SPC25","ST8SIA1","STAR","TBCA","TMCC3","TP53I11","VPS29","WNT10B","ZNF483")

unique(c(E2vsE3,E4vsE3,E4vsE2))

seurat_E2vsE3<-subset(x = seurat, features=E2vsE3)
seurat_E4vsE2<-subset(x = seurat, features=E4vsE2)
seurat_E4vsE3<-subset(x = seurat, features=E4vsE3)
seurat_all<-subset(x = seurat, features=c(E2vsE3,E4vsE3,E4vsE2))

E2vsE3_percent <- Matrix::colSums(seurat_E2vsE3@assays$RNA@counts)/Matrix::colSums(seurat@assays$RNA@counts)
E4vsE2_percent <- Matrix::colSums(seurat_E4vsE2@assays$RNA@counts)/Matrix::colSums(seurat@assays$RNA@counts)
E4vsE3_percent <- Matrix::colSums(seurat_E4vsE3@assays$RNA@counts)/Matrix::colSums(seurat@assays$RNA@counts)
All_percent <- Matrix::colSums(seurat_all@assays$RNA@counts)/Matrix::colSums(seurat@assays$RNA@counts)

seurat_E2vsE3<- AddMetaData(object = seurat_E2vsE3, metadata = E2vsE3_percent , col.name = "E2vsE3_percent")
seurat_E4vsE2<- AddMetaData(object = seurat_E4vsE2, metadata = E4vsE2_percent , col.name = "E4vsE2_percent")
seurat_E4vsE3<- AddMetaData(object = seurat_E4vsE3, metadata = E4vsE3_percent , col.name = "E4vsE3_percent")
seurat_all<- AddMetaData(object = seurat_all, metadata = All_percent , col.name = "All_percent")

Idents(seurat)<-"broad.cell.type"
levels(Idents(seurat))
table(Idents(seurat))  

#Dim Plot

DimPlot(seurat, reduction = "umap", label = TRUE)

# Feature Plot

FeaturePlot(seurat_all, features = "All_percent", cols = c("green", "red"), reduction = "umap", min.cutoff="q10",max.cutoff="q90")
FeaturePlot(seurat_E2vsE3, features = "E2vsE3_percent", cols = c("green", "red"), reduction = "umap", min.cutoff="q10",max.cutoff="q90")
FeaturePlot(seurat_E4vsE3, features = "E4vsE3_percent", cols = c("green", "red"), reduction = "umap", min.cutoff="q10",max.cutoff="q90")
FeaturePlot(seurat_E4vsE2, features = "E4vsE2_percent", cols = c("green", "red"), reduction = "umap", min.cutoff="q10",max.cutoff="q90")

# Find Markers

markers_Mic <- FindMarkers(seurat, ident.1 = "Mic", logfc.threshold = 0.25, min.pct = 0.1)
head(markers_Mic)
intersect(genes, rownames(markers_Mic))
markers_Mic[rownames(markers_Mic) %in% genes, ]

microglia_subset <- subset(seurat, idents = "Mic")
microglia_subset <- FindClusters(microglia_subset, resolution = 0.5)

microglia_subset$gene_signature <- rowMeans(FetchData(microglia_subset, vars = genes), na.rm = TRUE)

VlnPlot(microglia_subset, features = "gene_signature")

gene_expression <- rowMeans(FetchData(microglia_subset, vars = genes), na.rm = TRUE)

microglia_subset$gene_group <- ifelse(gene_expression > 0, "Expressed", "Not Expressed")

DimPlot(microglia_subset, reduction = "umap", group.by = "gene_group", cols = c("blue", "gray"))



dir.create("featureplot_mic", showWarnings = FALSE)

for (gene in genes) {
  p <- FeaturePlot(microglia_subset, features = gene, reduction = "umap") +
    labs(title = gene)  
  
  ggsave(filename = paste0("featureplot_mic/FeaturePlot_", gene, ".png"),
         plot = p, width = 5, height = 5, dpi = 300)
}

dir.create("vlnplot_mic", showWarnings = FALSE)

for (gene in genes) {
  p <- VlnPlot(microglia_subset, features = gene) +
    labs(title = gene) 
  
  ggsave(filename = paste0("vlnplot_mic/VlnPlot_", gene, ".png"),
         plot = p, width = 5, height = 5, dpi = 300)
}

# FindAllMarker

all_markers <- FindAllMarkers(seurat, logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE)
head(all_markers)
DimPlot(microglia_subset, reduction = "umap", label = TRUE)

# FindMarkers Astrocitos

markers_Ast <- FindMarkers(seurat, ident.1 = "Ast", logfc.threshold = 0.25, min.pct = 0.1)
head(markers_Ast)
ast_subset <- subset(seurat, idents = "Ast")
ast_subset <- FindClusters(ast_subset, resolution = 0.5)
DimPlot(ast_subset, reduction = "umap", label = TRUE)
markers_Ast[rownames(markers_Ast) %in% genes, ]

# FindMarkers Endoteliales
markers_End <- FindMarkers(seurat, ident.1 = "End", logfc.threshold = 0.25, min.pct = 0.1)
head(markers_End)


end_subset <- subset(seurat, idents = "End")
end_subset <- FindClusters(end_subset, resolution = 0.5)
DimPlot(end_subset, reduction = "umap", label = TRUE)

markers_End[rownames(markers_End) %in% genes, ]

# FindMarkers Neuronas Excitatorias

markers_Ex <- FindMarkers(seurat, ident.1 = "Ex", logfc.threshold = 0.25, min.pct = 0.1)
head(markers_Ex)
ex_subset <- subset(seurat, idents = "Ex")
ex_subset <- FindClusters(ex_subset, resolution = 0.5)
DimPlot(ex_subset, reduction = "umap", label = TRUE)
markers_Ex[rownames(markers_Ex) %in% genes, ]

# FindMarkers Neuronas Inhibitorias
markers_In <- FindMarkers(seurat, ident.1 = "In", logfc.threshold = 0.25, min.pct = 0.1)
head(markers_In)
in_subset <- subset(seurat, idents = "In")
in_subset <- FindClusters(in_subset, resolution = 0.5)
DimPlot(in_subset, reduction = "umap", label = TRUE)

markers_In[rownames(markers_In) %in% genes, ]

# FindMarkers Oligodendrocitos

markers_Oli <- FindMarkers(seurat, ident.1 = "Oli", logfc.threshold = 0.25, min.pct = 0.1)
head(markers_Oli)

oli_subset <- subset(seurat, idents = "Oli")
oli_subset <- FindClusters(oli_subset, resolution = 0.5)
DimPlot(oli_subset, reduction = "umap", label = TRUE)

markers_Oli[rownames(markers_Oli) %in% genes, ]

# FindMarkers Precursores Oligodendrocitos
markers_Opc <- FindMarkers(seurat, ident.1 = "Opc", logfc.threshold = 0.25, min.pct = 0.1)
head(markers_Opc)
opc_subset <- subset(seurat, idents = "Opc")
opc_subset <- FindClusters(opc_subset, resolution = 0.5)
DimPlot(opc_subset, reduction = "umap", label = TRUE)
markers_Opc[rownames(markers_Opc) %in% genes, ]

# FindMarkers Pericitos
markers_Per <- FindMarkers(seurat, ident.1 = "Per", logfc.threshold = 0.25, min.pct = 0.1)
head(markers_Per)
per_subset <- subset(seurat, idents = "Per")
per_subset <- FindClusters(per_subset, resolution = 0.5)
DimPlot(per_subset, reduction = "umap", label = TRUE)
markers_Per[rownames(markers_Per) %in% genes, ]

# For saving the markers in a excel

library(openxlsx)

wb <- createWorkbook()

markers_list <- list(
  "Ast" = markers_Ast[rownames(markers_Ast) %in% genes, ],
  "End" = markers_End[rownames(markers_End) %in% genes, ],
  "Ex" = markers_Ex[rownames(markers_Ex) %in% genes, ],
  "In" = markers_In[rownames(markers_In) %in% genes, ],
  "Mic" = markers_Mic[rownames(markers_Mic) %in% genes, ],
  "Oli" = markers_Oli[rownames(markers_Oli) %in% genes, ],
  "Opc" = markers_Opc[rownames(markers_Opc) %in% genes, ],
  "Per" = markers_Per[rownames(markers_Per) %in% genes, ]
)

for (name in names(markers_list)) {
  addWorksheet(wb, name)
  writeData(wb, sheet = name, markers_list[[name]], rowNames = TRUE)
}

saveWorkbook(wb, "Markers.xlsx", overwrite = TRUE)

# Monocle

seurat_matrix <- seurat[["RNA"]]$counts
cell_metadata <- seurat@meta.data
gene_metadata <- seurat[["RNA"]]@meta.features

library(monocle3)

cds <- new_cell_data_set(
  seurat_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

cds <- preprocess_cds(cds, method = "PCA")

cds <- reduce_dimension(cds, max_components = 2)

cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "cluster")


get_earliest_principal_node <- function(cds, time_bin="Mic"){
  cell_ids <- which(colData(cds)[, "broad.cell.type"] == time_bin)
  
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
plot_cells(cds, color_cells_by="broad.cell.type", group_cells_by="broad.cell.type", group_label_size = 10, label_groups = TRUE, show_trajectory_graph = TRUE, label_leaves = TRUE, label_branch_points=TRUE)


# Análisis por Genotipo

table(seurat@meta.data$ApoE_group)
DimPlot(seurat, group.by = "ApoE_group", reduction = "umap")
plot_cells(cds, color_cells_by = "ApoE_group")
markers_apoe_2_4 <- FindMarkers(seurat, ident.1 = "ApoE2", ident.2 = "ApoE4", group.by = "ApoE_group")
head(markers_apoe_2_4)
markers_apoe_3_4 <- FindMarkers(seurat, ident.1 = "ApoE3", ident.2 = "ApoE4", group.by = "ApoE_group")
head(markers_apoe_3_4)
markers_apoe_2_3 <- FindMarkers(seurat, ident.1 = "ApoE2", ident.2 = "ApoE3", group.by = "ApoE_group")

wb <- createWorkbook()

datasets <- list(
  "markers_apoe_2_4" = markers_apoe_2_4,
  "markers_apoe_3_4" = markers_apoe_3_4,
  "markers_apoe_2_3" = markers_apoe_2_3
)

for (name in names(datasets)) {
  addWorksheet(wb, name)
  writeData(wb, sheet = name, datasets[[name]])
}

saveWorkbook(wb, "markers_apoe.xlsx", overwrite = TRUE)


plot_cells(cds, color_cells_by = "pseudotime", 
           show_trajectory_graph = TRUE)

top_markers_res <- top_markers(cds, 
                               group_cells_by = "broad.cell.type", 
                               cores = 4)  

head(top_markers_res)


top_specific_markers <- top_markers_res %>% 
  filter(marker_test_q_value <= 0.05,  
         fraction_expressing >= 0.10) %>%  
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)  

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

top_specific_markers


rowData(cds)$gene_short_name <- rownames(rowData(cds))

plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by = "broad.cell.type",
                    ordering_type = "cluster_row_col",
                    max.size = 9)

plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by = "ApoE_group",
                    ordering_type = "cluster_row_col",
                    max.size = 9)
