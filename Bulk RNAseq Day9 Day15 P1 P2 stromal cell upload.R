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

mydata <- read.csv("P1P2counts.csv")
Gene_length <- data.frame(row.names = mydata[,6], mydata$Stop-mydata$Start)
colnames(Gene_length) <- "Length"
sampleinfo <- read.csv("sampleinfo.csv")

mydata_counts <- mydata[,-(1:9)]
rownames(mydata_counts) <- mydata[,6]
colnames(mydata_counts) <- sampleinfo$SampleName

mydata_TPM <- as.data.frame(matrix(0, 20252, 16))
for(j in 1:16){
  S <- 0
  for(i in 1:20252){
    S <- S + mydata_counts[i, j] / Gene_length[i, 1]
  }
  for(n in 1:20252){
    mydata_TPM[n,j] <- 10^6 * mydata_counts[n,j] / (as.numeric(Gene_length[n, 1]) * S)
  }
}
colnames(mydata_TPM) <-colnames(mydata_counts)
rownames(mydata_TPM) <-rownames(mydata_counts)
write.csv(mydata_TPM, file = "D9D15P1P2_TPM.csv")

keep_genes <- rowSums(mydata_counts == 0) <= 8
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
    x = "PC1 (21.1%)",
    y = "PC2 (13.6%)",
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
res_D9P1_vs_D15P1 <- results(dds, contrast = c("Status", "D9 P1", "D15 P1"))
res_D9P2_vs_D15P2 <- results(dds, contrast = c("Status", "D9 P2", "D15 P2"))
res_D9P1_vs_D9P2 <- results(dds, contrast = c("Status", "D9 P1", "D9 P2"))
summary(res_D9P1_vs_D15P1)
summary(res_D9P2_vs_D15P2)
summary(res_D9P1_vs_D9P2)

log2FC_threshold <- 1
padj_threshold <- 0.05

#groups
res <- res_D9P1_vs_D15P1
res <- res_D9P2_vs_D15P2
res <- res_D9P1_vs_D9P2

res <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]
summary(res)
DEGseg <- res[(abs(res$log2FoldChange) > log2FC_threshold) & (res$padj < padj_threshold), ]
DEGseg <- as.data.frame(DEGseg)
DEGseg$Gene <- rownames(DEGseg)
DEGseg_merged <- merge(DEGseg, mydata_TPM, by.x = "Gene", by.y = "row.names", all.x = TRUE)

plotCounts(dds, gene="Bmp2", intgroup="Status")
normalized_counts <- counts(dds, normalized = TRUE)
bmp2_counts <- normalized_counts["Bmp2", ]
bmp2_counts

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
genes_to_label <- c("Wnt4","Agt","Adamdec1","Ncam1","Dach1","Col6a4",
                    "Cd9","Cd34","Cmah","Igfbp6","Igf2","Efemp1")
label_genes <- subset(results_df, Gene %in% genes_to_label)

ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Color), alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", 
                                "Downregulated" = "blue", 
                                "Not Significant" = "gray")) +
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed", color = "black") +
  geom_text(data = label_genes, aes(label = Gene), size = 3, vjust = 1, hjust = 1) + #label genes
  xlim(-15, 15) +  
  theme_minimal() +
  labs(title = "D9 P1 vs D9 P2 Volcano Plot", 
       x = "Log2 Fold Change", 
       y = "-Log10 Adjusted P-Value") +
  theme(panel.background = element_rect(fill = "white"))

ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Color), alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", 
                                "Downregulated" = "blue", 
                                "Not Significant" = "gray")) +
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed", color = "black") +
  xlim(-12, 12) +  
  ylim(-0.5, 17) +  
  theme_minimal() +
  labs(title = "D9 P1 vs D15 P1 Volcano Plot", 
       x = "Log2 Fold Change", 
       y = "-Log10 Adjusted P-Value") +
  theme(panel.background = element_rect(fill = "white"))

ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Color), alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", 
                                "Downregulated" = "blue", 
                                "Not Significant" = "gray")) +
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed", color = "black") +
  xlim(-15, 15) +  
  ylim(-0.5, 17) + 
  theme_minimal() +
  labs(title = "D9 P2 vs D15 P2 Volcano Plot", 
       x = "Log2 Fold Change", 
       y = "-Log10 Adjusted P-Value") +
  theme(panel.background = element_rect(fill = "white"))

#heatmap
library(pheatmap)
selected_genes <- c("Wnt4","Agt","Syt13","Fam162b","Col15a1","Dach1","Adamdec1","Aard","Fhl2","Ctsc",
                    "Col6a4","Pdlim3","Fgd6","Ackr3","Tmem176a","Zfp536","Pitx1","Snai2","Rdh10","Rgs2")
present_genes <- selected_genes[selected_genes %in% rownames(mydata_TPM)]
selected_gene_data <- mydata_TPM[present_genes,]
norm_selected_gene_data <- log2(selected_gene_data + 1)
pheatmap(norm_selected_gene_data, ###FigS4b export 4.5*4.5
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         scale = "row",
         show_rownames = TRUE, 
         show_colnames = TRUE,
         main = "Heatmap of top 20 C1 Wnt4 stromal cell markers")

#intersect of DE gene list
library(VennDiagram)
DEP1 <- read.csv("DEGs_D9P1_vs_D15P1.csv", header = TRUE, stringsAsFactors = FALSE)
DEP2 <- read.csv("DEGs_D9P2_vs_D15P2.csv", header = TRUE, stringsAsFactors = FALSE)
genes_DEP1 <- DEP1[[1]]
genes_DEP2 <- DEP2[[1]]
unique_DEP1 <- setdiff(genes_DEP1, genes_DEP2)
venn.plot <- draw.pairwise.venn(
  area1 = length(genes_DEP1),
  area2 = length(genes_DEP2),
  cross.area = length(intersect(genes_DEP1, genes_DEP2)),
  category = c("DEP1", "DEP2"),
  fill = c("blue", "red"),
  alpha = 0.5,
  lty = "dashed",
  cex = 1.5,
  cat.cex = 1.2
)
unique_DEP1 <- as.data.frame(unique_DEP1)
unique_DEP1$Gene <- unique_DEP1[,1]
unique_DEP1_merged <- merge(unique_DEP1, DEP1, by.x = "Gene", by.y = "Gene", all.x = TRUE)
write.csv(unique_DEP1_merged, "Unique_DEP1_genes.csv", row.names = FALSE)
