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

mydata <- read.csv("wholetissue.csv")
Gene_length <- data.frame(row.names = mydata[,1], mydata$Stop-mydata$Start)
colnames(Gene_length) <- "Length"
sampleinfo <- read.csv("sampleinfo.csv")
mydata_counts <- mydata[,-(1:5)]
rownames(mydata_counts) <- mydata[,1]
colnames(mydata_counts) <- sampleinfo$SampleName

#PCA
keep_genes <- rowSums(mydata_counts == 0) <= 3
filtered_counts <- mydata_counts[keep_genes, ]
log_data <- log2(filtered_counts + 1)
group <- factor(sampleinfo$Status)
pca <- prcomp(t(log_data), scale. = FALSE)
summary(pca)
pca_df <- as.data.frame(pca$x)
pca_df$sample <- colnames(filtered_counts)
pca_df$group <- group

ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +  # Set the size of points
  labs(
    title = "PCA Plot",
    x = "PC1 (33.6%)",
    y = "PC2 (18.9%)",
    color = "Group"
  ) +
  theme_minimal(base_size = 10) +  # Minimal theme with larger text size
  theme(
    #panel.grid = element_blank(),  # Remove grid lines
    panel.background = element_rect(fill = "white"),  # White background
    #axis.line = element_line(color = "black"),  # Black axis lines
    plot.title = element_text(hjust = 0.5),  # Center the title
    legend.title = element_text(size = 10),  # Customize legend title
    legend.text = element_text(size = 10)    # Customize legend text
  )

#Heatmap
library(pheatmap)
selected_genes <- c("Pdgfra", "Lum", "Dcn", "Mfap5","Cxcl12","Col1a1","Col1a2","Col3a1","Col6a1","Col6a2","Col6a3","Agt","Adamdec1")
present_genes <- selected_genes[selected_genes %in% rownames(mydata_TPM)]
selected_gene_data <- mydata_TPM[present_genes,]
norm_selected_gene_data <- log2(selected_gene_data + 1)
pheatmap(norm_selected_gene_data, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         scale = "row",
         show_rownames = TRUE, 
         show_colnames = FALSE,
         main = "Heatmap of Stromal genes")


