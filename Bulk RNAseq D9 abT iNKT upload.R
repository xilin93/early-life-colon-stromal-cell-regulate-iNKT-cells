library(here)
library(edgeR)
library(Rtsne)
library(ggplot2)
library(plotly)
library(data.table)
library(DESeq2)
library(EnhancedVolcano)
library(RColorBrewer)
library(colorspace)
library(readr)

mydata <- read.csv("D9_NKT_abT.csv")
Gene_length <- data.frame(row.names = mydata[,1], mydata$Stop-mydata$Start)
colnames(Gene_length) <- "Length"
sampleinfo <- read.csv("sampleinfo.csv")
mydata_counts <- mydata[,-(1:8)]
rownames(mydata_counts) <- mydata[,1]
colnames(mydata_counts) <- sampleinfo$SampleName

mydata_TPM <- as.data.frame(matrix(0, 55291, 8))
for(j in 1:8){
  S <- 0
  for(i in 1:55291){
    S <- S + mydata_counts[i, j] / Gene_length[i, 1]
  }
  for(n in 1:55291){
    mydata_TPM[n,j] <- 10^6 * mydata_counts[n,j] / (as.numeric(Gene_length[n, 1]) * S)
  }
}
colnames(mydata_TPM) <-colnames(mydata_counts)
rownames(mydata_TPM) <-rownames(mydata_counts)

#PCA
keep_genes <- rowSums(mydata_counts == 0) <= 4
filtered_counts <- mydata_counts[keep_genes, ]
log_data <- log2(filtered_counts + 1)
group <- factor(sampleinfo$Status)
pca <- prcomp(t(log_data), scale. = FALSE)
summary(pca)
pca_df <- as.data.frame(pca$x)
pca_df$sample <- colnames(filtered_counts)
pca_df$group <- group

ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 2) +  
  stat_ellipse(
    aes(group = group, color = group),
    type = "norm",  
    linetype = "solid", 
    size = 1, 
    alpha = 0.3,  
    level = 0.6  
  ) +
  labs(
    title = "PCA Plot",
    x = "PC1 (27.2%)",
    y = "PC2 (21.5%)",
    color = "Group"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  ) 

#DEGs
is_integer <- all(filtered_counts == round(filtered_counts))
if (!is_integer) {
  filtered_counts <- round(filtered_counts)
}
sampleinfo$Status <- as.factor(sampleinfo$Status)
dds <- DESeqDataSetFromMatrix(countData = filtered_counts, 
                              colData = sampleinfo, 
                              design = ~ Status)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, contrast = c("Status", "D9 NKT", "D9 abT"))
summary(res)
log2FC_threshold <- 1
padj_threshold <- 0.05
res <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]
summary(res)
DEGseg <- res[(abs(res$log2FoldChange) > log2FC_threshold) & (res$padj < padj_threshold), ]
DEGseg <- as.data.frame(DEGseg)
DEGseg$Gene <- rownames(DEGseg)
DEGseg_merged <- merge(DEGseg, mydata_TPM, by.x = "Gene", by.y = "row.names", all.x = TRUE)
write.csv(DEGseg_merged, "DEGs_D9NKT_vs_D9abT.csv", row.names = FALSE)

#volcano plot 
results_df <- as.data.frame(res)
results_df$Significant <- with(results_df, 
                               ifelse(abs(log2FoldChange) > log2FC_threshold & padj < padj_threshold, 
                                      "Significant", 
                                      "Not Significant"))
results_df$Color <- with(results_df, 
                         ifelse(log2FoldChange > log2FC_threshold & padj < padj_threshold, 
                                "Upregulated", 
                                ifelse(log2FoldChange < -log2FC_threshold & padj < padj_threshold, 
                                       "Downregulated", 
                                       "Not Significant")))
results_df$Gene <- rownames(results_df)
label_genes <- subset(results_df, Significant == "Significant")
genes_to_label <- c("Ifitm1","Gzma","Mki67","Gzmb","Cxcr6","Tnf","Map2k3",
                    "Cd86","Foxp3","Lrrc32","Cd8b1","Ccr7","Il7r")
label_genes <- subset(results_df, Gene %in% genes_to_label)

ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Color), alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", 
                                "Downregulated" = "blue", 
                                "Not Significant" = "gray")) +
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed", color = "black") +
  geom_text(data = label_genes, aes(label = Gene), size = 3, vjust = 1, hjust = 1) + #label genes
  xlim(-15, 15) +  # Set the X-axis limits
  ylim(-0.5, 15) +  # Set the Y-axis limits
  theme_minimal() +
  labs(title = "D9 iNKT vs D9 abT Volcano Plot", 
       x = "Log2 Fold Change", 
       y = "-Log10 Adjusted P-Value") +
  theme(panel.background = element_rect(fill = "white"))

#Heatmap
library(pheatmap)
selected_genes <- c("Mapk6","Mapk7","Mapk9","Map2k3","Map2k4",
                    "Bmpr2","Acvr2a",
                    "Myc","Egr1","Rps6ka1","Cdt1","Dusp1","Dusp2","Lamtor3","Nfkb1","Nr4a1","Itgal","Icam1",
                    "Itga1","Itga2","Itga4","Itgb2","Itgb3","Itgb7")
present_genes <- selected_genes[selected_genes %in% rownames(mydata_TPM)]
selected_gene_data <- mydata_TPM[present_genes,]
norm_selected_gene_data <- log2(selected_gene_data + 1)
pheatmap(norm_selected_gene_data,
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         scale = "row",
         show_rownames = TRUE, 
         show_colnames = TRUE,
         main = "Heatmap")


