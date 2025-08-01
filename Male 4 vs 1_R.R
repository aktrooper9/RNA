# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(tidyr)
library("apeglm")
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ReactomePA)
library(umap)

# Create results directories
dir.create("results", showWarnings = FALSE)
dir.create("results/plots", showWarnings = FALSE)
dir.create("results/tables", showWarnings = FALSE)
dir.create("results/rds_plots", showWarnings = FALSE)


# Read processed metadata file
pdata <- read.csv("sample_metadata_gse115828.csv")

# Filter for males only and clean mgs_level
pdata_clean <- pdata[pdata$Sex == "M" & !pdata$mgs_level %in% c("no grade") & !is.na(pdata$mgs_level), ]

# Keep only samples with mgs_level "1" or "4"
pdata_clean <- pdata_clean[pdata_clean$mgs_level %in% c("1", "4"), ]
pdata_clean$mgs_level <- factor(pdata_clean$mgs_level, levels = c("1", "4"))

# Explicitly set mgs_level 1 as reference
pdata_clean$mgs_level <- relevel(pdata_clean$mgs_level, ref = "1")

# Check levels of mgs_level
print("mgs_level distribution in males (1 vs 4):")
print(levels(pdata_clean$mgs_level))
print(table(pdata_clean$mgs_level))

# Read counts file
counts_file <- "GSE115828_raw_counts_GRCh38.p13_NCBI.tsv"
counts <- read.table(counts_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# Check that sample IDs match between pdata_clean and counts column names
print(paste("All samples in metadata found in counts:", all(pdata_clean$sample_id %in% colnames(counts))))

# Match samples between metadata and counts
common_samples <- intersect(pdata_clean$sample_id, colnames(counts))
pdata_matched <- pdata_clean[pdata_clean$sample_id %in% common_samples, ]
pdata_matched <- pdata_matched[match(common_samples, pdata_matched$sample_id), ]

print(paste("Final sample count:", nrow(pdata_matched)))
print("Final mgs_level distribution:")
print(table(pdata_matched$mgs_level))

# Subset and reorder the counts matrix columns to match pdata_matched
counts_matched <- counts[, common_samples]

# Ensure rownames of metadata match colnames of counts
rownames(pdata_matched) <- pdata_matched$sample_id

# Create the DESeq2 dataset with mgs_level as the design
dds <- DESeqDataSetFromMatrix(countData = counts_matched,
                              colData = pdata_matched,
                              design = ~ mgs_level)

# Run DESeq2
dds <- DESeq(dds)

# Get results: mgs_level "4" vs "1" (4 is numerator, 1 is denominator)
res_mgs <- results(dds, contrast = c("mgs_level", "4", "1"))
norm_counts <- counts(dds, normalized = TRUE)

# Define thresholds
lfc_thresholds <- c(log2(1.5), log2(2))  # 0.58 and 1
fdr_threshold <- 0.05

# Function to filter and count DEGs for given thresholds
filter_deg_summary <- function(lfc_cutoff, fdr_cutoff, res) {
  filtered <- res[!is.na(res$padj) & 
                    res$padj < fdr_cutoff & 
                    abs(res$log2FoldChange) >= lfc_cutoff, ]
  n_up <- sum(filtered$log2FoldChange >= lfc_cutoff)
  n_down <- sum(filtered$log2FoldChange <= -lfc_cutoff)
  total <- nrow(filtered)
  
  data.frame(
    LFC_Cutoff = lfc_cutoff,
    FDR_Cutoff = fdr_cutoff,
    Upregulated = n_up,
    Downregulated = n_down,
    Total_DEGs = total
  )
}

# Apply for all LFC thresholds
results_summary <- do.call(rbind, lapply(lfc_thresholds, filter_deg_summary, fdr_cutoff = fdr_threshold, res = res_mgs))
write.csv(results_summary, "results/tables/deg_summary_mgs_1_vs_4_males.csv")

# Plot DEG counts
df_long <- pivot_longer(results_summary, cols = c("Upregulated", "Downregulated", "Total_DEGs"),
                        names_to = "Category", values_to = "Count")
df_long$LFC_Cutoff <- factor(df_long$LFC_Cutoff, labels = c("1.5 Fold", "2 Fold"))

png("results/plots/deg_counts_plot_mgs_1_vs_4_males.png", width = 800, height = 600)
ggplot(df_long, aes(x = LFC_Cutoff, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "DEG Counts: mgs_level 4 vs 1 (males Only)",
       x = "Fold Change Threshold", y = "Number of DEGs") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal()
dev.off()

# ------------------------
# Save DESeq2 Results
# ------------------------
write.csv(as.data.frame(res_mgs), "results/tables/deseq2_results_mgs_4_vs_1_males.csv")

res_sig <- res_mgs[!is.na(res_mgs$padj) & res_mgs$padj < fdr_threshold & abs(res_mgs$log2FoldChange) > log2(1.5), ]
write.csv(as.data.frame(res_sig), "results/tables/deseq2_sig_results_mgs_1_vs_4_males.csv")
write.csv(norm_counts, "results/tables/normalized_counts_mgs_1_vs_4_males.csv")

# ------------------------
# Plots: MA, Volcano, Top DEGs
# ------------------------
png("results/plots/ma_plot_mgs_1_vs_4_males.png", width = 800, height = 600)
plotMA(res_mgs, ylim = c(-4, 4), main = "MA Plot: mgs_level 4 vs 1 (males)")
dev.off()

res_df <- as.data.frame(res_mgs)
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > log2(1.5), "Yes", "No")

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  labs(title = "Volcano Plot: mgs_level 4 vs 1 (males)",
       x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme_minimal()

ggsave("results/plots/volcano_plot_mgs_1_vs_4_males.png", volcano_plot, width = 8, height = 6)
saveRDS(volcano_plot, file = "results/rds_plots/volcano_plot_mgs_1_vs_4_males.rds")

# Top 10 downregulated (lowest log2FoldChange)
top_down <- head(res_sig[order(res_sig$log2FoldChange), ], 10)

# Top 10 upregulated (highest log2FoldChange)
top_up <- head(res_sig[order(-res_sig$log2FoldChange), ], 10)

# Combine them into one data frame
top_genes <- rbind(top_down, top_up)
if(nrow(top_genes) > 0) {
  top_degs_plot <- ggplot(top_genes, aes(x = reorder(rownames(top_genes), log2FoldChange),
                                         y = log2FoldChange, fill = log2FoldChange > 0)) +
    geom_col(show.legend = FALSE) +
    coord_flip() +
    labs(title = "Top 10 DEGs: mgs_level 4 vs 1 (males)",
         x = "Gene", y = "Log2 Fold Change") +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
    theme_minimal()
  ggsave("results/plots/top_degs_plot_mgs_1_vs_4_males.png", top_degs_plot, width = 8, height = 6)
  saveRDS(top_degs_plot, file = "results/rds_plots/top_degs_plot_mgs_1_vs_4_males.rds")
}

# ------------------------
# Enrichment Analysis
# ------------------------

sig_entrez <- rownames(res_sig)
bg_entrez <- rownames(res_mgs)

ego <- enrichGO(gene = sig_entrez, universe = bg_entrez,
                OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

ekegg <- enrichKEGG(gene = sig_entrez, universe = bg_entrez,
                    organism = "hsa", pvalueCutoff = 0.05)

react <- enrichPathway(gene = sig_entrez, universe = bg_entrez,
                       organism = "human", pvalueCutoff = 0.05, readable = TRUE)

# GO Plot
go_dot <- dotplot(ego, showCategory = 20) + ggtitle("GO: Biological Process Enrichment")
png("results/plots/GO_BP_enrichment_1_vs_4.png", width = 1000, height = 800)
print(go_dot)
dev.off()
saveRDS(go_dot, file = "results/rds_plots/go_bp_dotplot_mgs_1_vs_4_males.rds")

# KEGG Plot
kegg_dot <- dotplot(ekegg, showCategory = 20) + ggtitle("KEGG Pathway Enrichment")
png("results/plots/KEGG_enrichment_1_vs_4.png", width = 1000, height = 800)
print(kegg_dot)
dev.off()
saveRDS(kegg_dot, file = "results/rds_plots/kegg_dotplot_mgs_1_vs_4_males.rds")

# Reactome Plot
react_bar <- barplot(react, showCategory = 20) + ggtitle("Reactome Pathway Enrichment")
png("results/plots/Reactome_enrichment_1_vs_4.png", width = 1000, height = 800)
print(react_bar)
dev.off()
saveRDS(react_bar, file = "results/rds_plots/reactome_barplot_mgs_1_vs_4_males.rds")

write.csv(as.data.frame(ego), "results/tables/go_bp_enrichment_1_vs_4.csv")
write.csv(as.data.frame(ekegg), "results/tables/kegg_enrichment_1_vs_4.csv")
write.csv(as.data.frame(react), "results/tables/reactome_enrichment_1_vs_4.csv")

# ------------------------
# Heatmap of Top DEGs
# ------------------------

top_genes_heatmap <- if(nrow(res_sig) >= 30) rownames(head(res_sig[order(res_sig$padj), ], 30)) else rownames(res_sig)

if(length(top_genes_heatmap) > 0) {
  top_counts <- norm_counts[top_genes_heatmap, ]
  top_counts_scaled <- t(scale(t(top_counts)))
  
  annotation_col <- as.data.frame(colData(dds)[, "mgs_level", drop = FALSE])
  sample_order <- rownames(annotation_col)[order(annotation_col$mgs_level)]
  top_counts_sorted <- top_counts_scaled[, sample_order]
  annotation_col_sorted <- annotation_col[sample_order, , drop = FALSE]
  
  pheatmap(top_counts_sorted,
           cluster_rows = TRUE, cluster_cols = FALSE,
           annotation_col = annotation_col_sorted,
           show_rownames = TRUE, show_colnames = TRUE,
           fontsize_row = 8, fontsize_col = 9,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
           main = paste("Heatmap of Top", length(top_genes_heatmap), "DEGs: mgs_level 4 vs 1 (males)"),
           filename = "results/plots/heatmap_top_DEGs_mgs_1_vs_4_males.png")
  
  saveRDS(top_counts_sorted, "results/rds_plots/heatmap_matrix_top_DEGs_mgs_1_vs_4.rds")
}

# ------------------------
# UMAP
# ------------------------
variances <- apply(norm_counts, 1, var)
filtered_normalized_counts <- norm_counts[variances > 0, ]
normalized_counts_zscore <- scale(t(filtered_normalized_counts))
umap_result <- umap(normalized_counts_zscore)

umap_data <- as.data.frame(umap_result$layout)
umap_data$sample_id <- rownames(umap_data)
umap_data$mgs_level <- pdata_matched$mgs_level[match(umap_data$sample_id, pdata_matched$sample_id)]

umap_plot <- ggplot(umap_data, aes(x = V1, y = V2, color = mgs_level)) +
  geom_point(size = 3) +
  labs(title = "UMAP of Normalized Counts: mgs_level Groups (males Only)",
       x = "UMAP1", y = "UMAP2") +
  theme_minimal() +
  scale_color_manual(values = c("1" = "blue", "4" = "red"))

ggsave("results/plots/umap_plot_mgs_1_vs_4_males.png", umap_plot, width = 8, height = 6)
saveRDS(umap_plot, file = "results/rds_plots/umap_plot_mgs_1_vs_4_males.rds")

# ------------------------
# Session Info and Final Messages
# ------------------------
writeLines(capture.output(sessionInfo()), "results/session_info_mgs_1_vs_4_males.txt")

cat("\nAnalysis complete!\n")
cat("Total samples analyzed:", nrow(pdata_matched), "\n")
cat("mgs_level group sizes:\n")
print(table(pdata_matched$mgs_level))

cat("\nSummary of significant DESeq2 results (res_sig):\n")
print(summary(res_sig))
