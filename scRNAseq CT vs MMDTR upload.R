
library(Seurat)
library(SeuratObject)
library(dplyr)
library(here)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(RColorBrewer)
library(viridis)

# load data
wt_data <- Read10X(data.dir = here("WT_filtered_matrix"))
mmdtr_data <- Read10X(data.dir = here("MMDTR_filtered_matrix"))
CT <- CreateSeuratObject(counts = wt_data, project = "CT", min.cells = 3, min.features = 200)
MMDTR <- CreateSeuratObject(counts = mmdtr_data, project = "MMDTR", min.cells = 3, min.features = 200)
CT$group <- "CT"
MMDTR$group <- "MMDTR"

# Quality Control
CT[["percent.mt"]] <- PercentageFeatureSet(CT, pattern = "^mt-")
MMDTR[["percent.mt"]] <- PercentageFeatureSet(MMDTR, pattern = "^mt-")
VlnPlot(CT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(MMDTR, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CT <- subset(CT, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
MMDTR <- subset(MMDTR, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)

# Normalize and Data Integration
CT <- NormalizeData(CT)
CT <- FindVariableFeatures(CT, selection.method = "vst", nfeatures = 2000)
MMDTR <- NormalizeData(MMDTR)
MMDTR <- FindVariableFeatures(MMDTR, selection.method = "vst", nfeatures = 2000)
anchors <- FindIntegrationAnchors(object.list = list(CT, MMDTR), dims = 1:20)
combined <- IntegrateData(anchorset = anchors)

# Scale Data Clustering and UMAP
combined <- ScaleData(combined)
combined <- RunPCA(combined, npcs = 30)
combined <- FindNeighbors(combined, dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)
combined <- RunUMAP(combined, dims = 1:20)

# Step 7: Identify Cell Type Markers
cluster_markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster_markers, "cluster_markers.csv")

saveRDS(combined, file = here("combined.rds"))

# rename each cluster
new_cluster_names <- c(
  "0" = "C0 Angptl1 stromal cell",
  "1" = "C1 Wnt4 stromal cell",
  "2" = "C2 Fap stromal cell",
  "3" = "C3 Adgrd1 stromal cell",
  "4" = "C4 Proliferating cell",
  "5" = "C5 Igf2 stromal cell",
  "6" = "C6 Chodl stromal cell",
  "7" = "C7 Glia cell",
  "8" = "C8 Fgf9 stromal cell",
  "9" = "C9 ILC",
  "10" = "C10 Dendritic Cell",
  "11" = "C11 Smooth muscle cell",
  "12" = "C12 Pericyte",
  "13" = "C13 Macrophage/Neutrophil",
  "14" = "C14 T cell",
  "15" = "C15 Endothelial cell")
rename <- RenameIdents(combined, new_cluster_names)
table(Idents(rename))

saveRDS(rename, file = here("rename.rds"))

DimPlot(rename, reduction = "umap", split.by="group", label = TRUE, label.size =5) + ggtitle("UMAP - CT vs MMDTR")
vln_plot <- VlnPlot(rename, features = c("Wnt4", "Agt", "Adamdec1"), pt.size = 0, ncol = 1)
vln_plot <- vln_plot + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate X-axis labels by 45 degrees
vln_plot#export 11*12 inch
FeaturePlot(combined, features = c("Wnt4","Pdgfra","Ncam1","Cd34","Cd9"), label = FALSE, ncol = 5)#export 3*16 inch
FeaturePlot(combined, features = c("Bmp2"), label = FALSE)#export 3*3.5 inch
FeaturePlot(combined, features = c("Foxl1"), label = FALSE)#export 3*3.5 inch
FeaturePlot(combined, features = c("Cd1d1"), label = FALSE)#export 3*3.5 inch

# subset Mac/neutro C13 cluster to analyze
macneu <- subset(combined, idents = "13")
table(macneu@meta.data$group)
macneu <- NormalizeData(macneu)
macneu <- FindVariableFeatures(macneu, selection.method = "vst", nfeatures = 2000)
macneu <- ScaleData(macneu)
macneu <- RunPCA(macneu, npcs = 30)
macneu <- FindNeighbors(macneu, dims = 1:20)
macneu <- FindClusters(macneu, resolution = 0.1)
macneu <- RunUMAP(macneu, dims = 1:20)
table(macneu@active.ident, macneu@meta.data$group)
macneu <- JoinLayers(macneu)
macneu_markers <- FindAllMarkers(macneu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(macneu_markers, "macneu_markers.csv")
new_cluster_names <- c("0" = "Macrophage","1" = "Neutrophil")
renamemacneu <- RenameIdents(macneu, new_cluster_names)
table(Idents(renamemacneu))
table(renamemacneu@active.ident, macneu@meta.data$group)

saveRDS(renamemacneu, file = here("Macneu.rds"))

DimPlot(renamemacneu, reduction = "umap", group.by = "group", label = TRUE, label.size =0) + ggtitle("UMAP_Macrophage/Neutrophil") # export 4*5inch
VlnPlot(renamemacneu, features = c("Adgre1", "Fcgr1", "Cx3cr1", "Csf3r","S100a8","S100a9"), group.by = "seurat_clusters", pt.size = 0.5)+xlab("")# export 4.5*6inch

# Differential Expression Analysis Between whole Groups
DefaultAssay(combined) <- "RNA"
combined <- JoinLayers(combined)
degs <- FindMarkers(combined, ident.1 = "CT", ident.2 = "MMDTR", group.by = "group")
write.csv(degs, "differential_genes_CT_vs_MMDTR.csv")
degs$significant <- ifelse(
  degs$p_val_adj < 0.05 & abs(degs$avg_log2FC) > 1,
  ifelse(degs$avg_log2FC > 1, "Downregulated in MMDTR", "Upregulated in MMDTR"),
  "Not Significant")
top_upregulated <- degs[degs$significant == "Upregulated in MMDTR", ]
top_downregulated <- degs[degs$significant == "Downregulated in MMDTR", ]
top_genes_to_label <- rbind(
  head(top_upregulated[order(-top_upregulated$avg_log2FC), ], 20),  # Top 20 upregulated
  head(top_downregulated[order(top_downregulated$avg_log2FC), ], 20)  # Top 20 downregulated
)
genes_to_label <- rownames(top_genes_to_label)
genes_to_label <- c("C1qa", "C1qb", "Cx3cr1","Fcgr1","Mrc1","Ms4a7","Fcrls",
                    "Stfa2l1","S100a8","S100a9","Ccl2","Ccl7","Cxcl2","Cxcl5","C3","Trem1")
degs$gene_label <- ifelse(rownames(degs) %in% genes_to_label, rownames(degs), "")
ggplot(degs, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significant)) +
  geom_point(alpha = 0.7) +
  geom_text(aes(label = gene_label), vjust =1.5, hjust = 0.8, size = 4, color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  scale_color_manual(
    values = c("blue","gray","red"),
    labels = c("Downregulated in MMDTR", "Not Significant", "Upregulated in MMDTR")
  ) +
  scale_y_continuous(limits = c(0, 20)) +
  scale_x_continuous(limits = c(-6, 6)) +
  theme_minimal(base_size = 15) +  # Adjust font size for better readability
  theme(
    panel.background = element_rect(fill = "white", color = "black"),  # White background
    panel.grid.major = element_line(color = "gray90"),  # Light grid lines
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Volcano Plot: CT vs MMDTR",
    x = "Log2 Fold Change (CT vs MMDTR)",
    y = "-Log10 Adjusted P-Value",
    color = "Gene Regulation")

feature_plots <- FeaturePlot(combined,features = c("Ptprc", "Coro1a", "Pecam1", "Cdh5", "Pdgfrb", "Cox4i2", 
                                                   "Foxd3", "Ptprz1", "Pdgfra", "Ptch1", "Lum", "Col6a1"), label = TRUE, label.size =3, ncol = 6)
feature_plots

dot_plot <- DotPlot(combined, features = c("Angptl1","Hmcn2","Grem1","C1qtnf3",
                                           "Wnt4","Agt","Adamdec1","Col15a1",
                                           "Fap","Ly6h","Shisa3","Mmp16",
                                           "Adgrd1","Cmah","Igfbp6","Efemp1",
                                           "Mki67","Cdk1","Ccnb1","Cdca8",
                                           "Igf2","Meox2","Prss35","Gata6",
                                           "Chodl","Ackr4","Fxyd6","Kcnn3",
                                           "Foxd3","Ptprz1","Nell2","Gpr37l1",
                                           "Fgf9","Dsp","Inhba","Bmp7",
                                           "Il7r","Gpr183","Il22","Il17f","Rorc","Klrb1b",
                                           "Cd209a","Flt3","Itgax","H2-Aa","Cd86",
                                           "Actg2","Cnn1","Lmod1","Des",
                                           "Map3k7cl","Notch3","Pdgfrb","Rgs5","Cox4i2",
                                           "Fcgr1","Cx3cr1","Csf1r","Csf3r",
                                           "Cd3g","Trbc1","Il2ra","Cxcr6",
                                           "Pecam1","Lyve1","Cdh5","Myct1")) +
  scale_color_gradient(low = "gray", high = "red") +  # Set color gradient
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate X-axis labels
  labs(x = "Genes", y = "Cluster Identity", color = "Expression") +  # Customize axis labels
  coord_flip() # Flip the axes
dot_plot