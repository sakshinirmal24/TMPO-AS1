# Load required libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)

# Step 1: Load the raw read count data (Replace 'your_file.csv' with your actual file)
counts <- read.csv("your_file.csv", row.names = 1, header = TRUE)

# Step 2: Load sample metadata (Ensure sample IDs match column names in counts matrix)
metadata <- read.csv("metadata.csv", row.names = 1, header = TRUE)

# Step 3: Filter only normal and LUAD samples
metadata <- metadata %>% filter(sample_type %in% c("Normal", "LUAD"))
counts <- counts[, colnames(counts) %in% rownames(metadata)]

# Step 4: Ensure column names of counts match metadata row names
metadata <- metadata[colnames(counts),]

# Step 5: Convert sample type to factor (DESeq2 requires factor format)
metadata$sample_type <- factor(metadata$sample_type, levels = c("Normal", "LUAD"))

# Step 6: Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ sample_type)

# Step 7: Normalize counts and run DESeq2 analysis
dds <- DESeq(dds)
res <- results(dds, contrast = c("sample_type", "LUAD", "Normal"))  # LUAD vs. Normal comparison

# Step 8: Filter for hsa-let-7b expression
let7b_results <- res[rownames(res) == "hsa-let-7b", ]
print(let7b_results)  # View results

# Step 9: Save results
write.csv(as.data.frame(res), file = "DESeq2_results.csv")

# Step 10: Generate Volcano Plot
res_df <- as.data.frame(res)
res_df$Significance <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Significant", "Not Significant")

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("red", "black")) +
  theme_minimal() +
  ggtitle("Volcano Plot: LUAD vs Normal") +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted p-value")

# Step 11: Boxplot for hsa-let-7b expression
norm_counts <- counts(dds, normalized = TRUE)
let7b_data <- data.frame(Expression = norm_counts["hsa-let-7b", ], Sample_Type = metadata$sample_type)

ggplot(let7b_data, aes(x = Sample_Type, y = Expression, fill = Sample_Type)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("hsa-let-7b Expression in LUAD vs Normal") +
  xlab("Sample Type") +
  ylab("Normalized Expression")

# Step 12: Heatmap for visualization
top_miRNAs <- rownames(res_df[order(res_df$padj), ])[1:20]  # Select top 20 differentially expressed miRNAs
pheatmap(norm_counts[top_miRNAs, ], cluster_rows = TRUE, cluster_cols = TRUE, scale = "row",
         annotation_col = metadata, show_rownames = TRUE, show_colnames = FALSE)
