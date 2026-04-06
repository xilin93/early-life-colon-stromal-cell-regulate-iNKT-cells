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

mydata <- read.csv("MMDTRNKT_count.csv")
Gene_length <- data.frame(row.names = mydata[,6], mydata$Stop-mydata$Start)
colnames(Gene_length) <- "Length"
sampleinfo <- read.csv("sampleinfo.csv")
mydata_counts <- mydata[,-(1:9)]
rownames(mydata_counts) <- mydata[,6]
colnames(mydata_counts) <- sampleinfo$SampleName

mydata_TPM <- as.data.frame(matrix(0, 16367, 8))
for(j in 1:8){
  S <- 0
  for(i in 1:16367){
    S <- S + mydata_counts[i, j] / Gene_length[i, 1]
  }
  for(n in 1:16367){
    mydata_TPM[n,j] <- 10^6 * mydata_counts[n,j] / (as.numeric(Gene_length[n, 1]) * S)
  }
}
colnames(mydata_TPM) <-colnames(mydata_counts)
rownames(mydata_TPM) <-rownames(mydata_counts)
write.csv(mydata_TPM, file = "CTMMDTR_NKT_TPM.csv")

#DEGs
keep_genes <- rowSums(mydata_counts == 0) <= 4
filtered_counts <- mydata_counts[keep_genes, ]
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
res <- results(dds, contrast = c("Status", "CT NKT", "MMDTR NKT"))
summary(res)
log2FC_threshold <- 1
padj_threshold <- 0.05
res <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]
summary(res)
DEGseg <- res[(abs(res$log2FoldChange) > log2FC_threshold) & (res$padj < padj_threshold), ]
DEGseg <- as.data.frame(DEGseg)
DEGseg$Gene <- rownames(DEGseg)
DEGseg_merged <- merge(DEGseg, mydata_TPM, by.x = "Gene", by.y = "row.names", all.x = TRUE)
write.csv(DEGseg_merged, "DEGs_CTNKT_vs_MMDTRNKT.csv", row.names = FALSE)

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
genes_to_label <- c("Cdc25b","Cit","Dusp6","Sh2b3",
                    "Rnf6","Ikzf5","Map7d1","Coa5")###FigS4k
label_genes <- subset(results_df, Gene %in% genes_to_label)

ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Color), alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", 
                                "Downregulated" = "blue", 
                                "Not Significant" = "gray")) +
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed", color = "black") +
  geom_text(data = label_genes, aes(label = Gene), size = 3, vjust = 1, hjust = 1) + #label genes
  xlim(-20, 20) +  # Set the X-axis limits
  ylim(-0.5, 12) +  # Set the Y-axis limits
  theme_minimal() +
  labs(title = "CT iNKT vs MMDTR iNKT Volcano Plot", 
       x = "Log2 Fold Change", 
       y = "-Log10 Adjusted P-Value") +
  theme(panel.background = element_rect(fill = "white"))

#Heatmap
library(pheatmap)
selected_genes <- c("Cdc25b","Gadd45a","Dusp6","Atf3","Fam162a","Trim39","Ccar2","Runx2","Fbxo5",
                    "Kps6ka1","Mapk14","Map2k1","Map2k4","Map2k5","Map3k2","Map3k7")
present_genes <- selected_genes[selected_genes %in% rownames(mydata_TPM)]
selected_gene_data <- mydata_TPM[present_genes,]
norm_selected_gene_data <- log2(selected_gene_data + 1)
pheatmap(norm_selected_gene_data,
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         scale = "row",
         show_rownames = TRUE, 
         show_colnames = TRUE,
         main = "Heatmap of top 20 C1 Wnt4 stromal cell markers")
