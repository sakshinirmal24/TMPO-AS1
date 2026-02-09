BiocManager::install(c("TCGAbiolinks", "tidyverse", "maftools", "pheatmap", "dplyr", "biomaRt", "ggpubr"))

library(TCGAbiolinks)
library(tidyverse)
library(maftools)
library(pheatmap)
library(SummarizedExperiment)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(ggpubr)

# get a list of projects
gdcprojects <- getGDCprojects()
getProjectSummary('TCGA-LUAD')

# building a query
query_TCGA <- GDCquery(project = 'TCGA-LUAD',
                       data.category = 'Transcriptome Profiling')
output_query_TCGA <- getResults(query_TCGA)

# build a query to retrieve gene expression data ------------
query_TCGA <- GDCquery(project = 'TCGA-LUAD',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       access = 'open')
getResults(query_TCGA)

# download data - GDCdownload
GDCdownload(query_TCGA)

# prepare data
tcga_luad_data <- GDCprepare(query_TCGA, summarizedExperiment = TRUE)
luad_matrix <- assay(tcga_luad_data, 'fpkm_unstrand')

# Extract clinical data
clinical_data <- as.data.frame(colData(tcga_luad_data)) 

# Select relevant columns
clinical_info <- clinical_data %>% dplyr::select(sample_type)

# Add sample IDs
clinical_info$sample_id <- rownames(clinical_info)

# View unique sample types
unique(clinical_info$sample_type)

# Extract STAR - Counts expression matrix
expr_matrix <- as.data.frame(assay(tcga_luad_data, "unstranded"))

# Convert rownames to column for gene names
expr_matrix <- rownames_to_column(expr_matrix, var = "gene")

# To identify TMPO-AS1 Ensembl ID
# Load Ensembl dataset
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get Ensembl ID for TMPO-AS1
TMPOS1_info <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = "TMPO-AS1",
  mart = mart
)

print(TMPOS1_info)

#Extract TMPO-AS1 Expression 
TMPOS1_expr <- expr_matrix[rownames(expr_matrix) == "ENSG00000257167", ]

#Remove Version Numbers
expr_matrix$gene <- sub("\\..*", "", expr_matrix$gene)
TMPOS1_expr <- expr_matrix %>% filter(gene == "ENSG00000257167")

#Extract TMPO-AS1 Expression
TMPOS1_expr_long <- as.data.frame(t(TMPOS1_expr[-1]))  # Remove gene column
colnames(TMPOS1_expr_long) <- "TMPOS1_expression"
TMPOS1_expr_long$sample_id <- rownames(TMPOS1_expr_long)

#Transpose for Merging
clinical_data$sample_id <- rownames(clinical_data)  # Ensure matching sample IDs
final_data <- merge(TMPOS1_expr_long, clinical_data, by = "sample_id")

filtered_data <- final_data %>%
  filter(sample_type %in% c("Primary Tumor", "Solid Tissue Normal"))

set.seed(123)  # For reproducibility
primary_subset <- filtered_data %>%
  filter(sample_type == "Primary Tumor") %>%
  sample_n(59)

normal_subset <- filtered_data %>%
  filter(sample_type == "Solid Tissue Normal")

colnames(balanced_data)[colnames(balanced_data) == "TMPO-AS1_expression"] <- "TMPOS1_expression"

balanced_data <- bind_rows(primary_subset, normal_subset)

table(balanced_data$sample_type)

stat_test <- wilcox.test(TMPOS1_expression ~ sample_type, data = balanced_data)
p_value <- stat_test$p.value  # Extract p-value

# Format p-value for display
p_value_text <- paste("p =", formatC(p_value, format = "e", digits = 2))


ggplot(balanced_data, aes(x = sample_type, y = TMPOS1_expression, fill = sample_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.5) +
  scale_fill_manual(values = c("blue", "red")) +  # Adjust colors as needed
  labs(title = "Expression of TMPO-AS1 in LUAD based on Sample Types",
       x = "TCGA Samples",
       y = "Transcript per million") +
  theme_minimal() +
  theme(text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(face = "bold"),
        legend.position = "none") + # Remove Legend
  
  # Add sample numbers below each category
  annotate("text", x = 1, y = min(balanced_data$TMPOS1_expression) - 200, 
           label = "(n=59)", size = 5, fontface = "bold") +
  annotate("text", x = 2, y = min(balanced_data$TMPOS1_expression) - 200, 
           label = "(n=59)", size = 5, fontface = "bold") +
  
  # Add p-value above the plot
  annotate("text", x = 1.5, y = max(balanced_data$TMPOS1_expression) * 1.1, 
           label = p_value_text, size = 6, fontface = "bold") +
  
  # Ensure equal spacing
  scale_x_discrete(labels = c("Primary Tumor", "Normal Tissue"))
