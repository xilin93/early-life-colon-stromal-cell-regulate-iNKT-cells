
library(Seurat)
library(SeuratObject)
library(Matrix)
library(dplyr)
library(here)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(RColorBrewer)
library(viridis)

barcodes_path <- here("GSE172261", "barcodes.tsv")
genes_path <- here("GSE172261", "genes.tsv")
matrix_path <- here("GSE172261", "matrix.mtx")
barcodes <- read.table(barcodes_path, stringsAsFactors = FALSE, header = FALSE)
genes <- read.table(genes_path, stringsAsFactors = FALSE, header = FALSE)
matrix <- readMM(matrix_path)
prefixes <- c("M1_Distal_H2O", "M1_Proximal_H2O", "M2_Distal_H2O", 
              "M2_Proximal_H2O", "M3_Distal_H2O", "M3_Proximal_H2O")
selected_barcodes <- barcodes %>%
  filter(sapply(V1, function(x) any(startsWith(x, prefixes))))
selected_indices <- which(barcodes$V1 %in% selected_barcodes$V1)
subset_matrix <- matrix[, selected_indices]
subset_barcodes <- barcodes[selected_indices, ]
write.table(subset_barcodes, "subset_barcodes.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
writeMM(subset_matrix, "subset_matrix.mtx")
write.table(genes, "subset_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

Day9_data <- Read10X(data.dir = here("WT_filtered_matrix"))
Adult_data <- Read10X(data.dir = here("GSE172261","GSE172261_subset"))
Day9 <- CreateSeuratObject(counts = Day9_data, project = "Day9", min.cells = 3, min.features = 200)
Adult <- CreateSeuratObject(counts = Adult_data, project = "Adult", min.cells = 3, min.features = 200)
Day9$group <- "Day9"
Adult$group <- "Adult"
Day9[["percent.mt"]] <- PercentageFeatureSet(Day9, pattern = "^mt-")
Adult[["percent.mt"]] <- PercentageFeatureSet(Adult, pattern = "^mt-")
Day9 <- subset(Day9, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
Adult <- subset(Adult, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 5)
Adult <- subset(Adult, downsample = 2000)

Day9 <- NormalizeData(Day9)
Day9 <- FindVariableFeatures(Day9, selection.method = "vst", nfeatures = 2000)
Adult <- NormalizeData(Adult)
Adult <- FindVariableFeatures(Adult, selection.method = "vst", nfeatures = 2000)
anchors <- FindIntegrationAnchors(object.list = list(Day9, Adult), dims = 1:20) #take time
combined <- IntegrateData(anchorset = anchors)

combined <- ScaleData(combined)
combined <- RunPCA(combined, npcs = 30)
combined <- FindNeighbors(combined, dims = 1:20)
combined <- FindClusters(combined, resolution = 0.3) #0.5-17, 0.3-15, 
combined <- RunUMAP(combined, dims = 1:20)

cluster_markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster_markers, "cluster_markers.csv")

saveRDS(combined, file = here("AdultDay9.rds"))

DimPlot(combined, reduction = "umap", group.by = "group", label = TRUE,label.size =0) + ggtitle("UMAP - Adult vs Day 9")
DimPlot(combined, reduction = "umap", split.by="group", label = TRUE,label.size =5) + ggtitle("UMAP - Adult vs Day 9")
DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE,label.size =5) + ggtitle("UMAP Clusters")
vln_plot <- VlnPlot(combined, features = c("Wnt4", "Agt"), pt.size = 0, ncol = 1)
vln_plot <- vln_plot + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
vln_plot
FeaturePlot(combined, features = c("Wnt4","Agt"), label = F, ncol = 2)h
FeaturePlot(combined, features = c("Wnt4", "Agt"), blend = TRUE) +
  ggtitle("Co-Expression of Wnt4 and Agt")
feature_plots <- FeaturePlot(combined,features = c("Pdgfra", "Ptch1", "Lum", "Dcn","Col6a1","Mki67",
                                                   "Actg2","Pecam1","Gpr37l1","Ptprc","Cspg4","Ano1"), label = TRUE, label.size =3, ncol = 6)  # Specify 6 columns
feature_plots

vln_plot <- VlnPlot(combined, features = c("Pdgfra", "Actg2","Ptch1", "Pecam1","Lum", "Gpr37l1","Dcn","Ptprc","Col6a1","Cspg4","Mki67","Ano1"), pt.size = 0, ncol = 2)
vln_plot <- vln_plot + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate X-axis labels
  )
vln_plot

clusters_of_interest <- c(0, 1, 2, 3, 4, 7, 8)
metadata <- combined@meta.data
wnt4_expression <- FetchData(combined, vars = "Wnt4")
metadata$Wnt4 <- wnt4_expression$Wnt4
subset_metadata <- metadata[metadata$seurat_clusters %in% clusters_of_interest & 
                              metadata$Wnt4 > 2, ]
table_by_group <- table(subset_metadata$group, subset_metadata$seurat_clusters)
print(table_by_group)


