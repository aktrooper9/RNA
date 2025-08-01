# ------------------------
# Setup and Libraries
# ------------------------
library(DESeq2)
library(ggplot2)
library(tidyr)
library(apeglm)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ReactomePA)
library(umap)
library(AnnotationDbi)

# Create results directories
dir.create("results_all", showWarnings = FALSE)
dir.create("results_all/plots", showWarnings = FALSE)
dir.create("results_all/tables", showWarnings = FALSE)
dir.create("results_all/rds_plots", showWarnings = FALSE)

# ------------------------
# Load Metadata and Filter (no sex filter)
# ------------------------
pdata <- read.csv("sample_metadata_gse115828.csv")
pdata_clean <- subset(pdata, mgs_level %in% c("1", "4") & !is.na(mgs_level))
pdata_clean$mgs_level <- factor(pdata_clean$mgs_level, levels = c("1", "4"))

# ------------------------
# Load and Align Counts
# ------------------------



common_samples <- intersect(pdata_clean$sample_id, colnames(counts))
pdata_matched <- pdata_clean[match(common_samples, pdata_clean$sample_id), ]
counts_matched <- counts[, common_samples]
rownames(pdata_matched) <- pdata_matched$sample_id

# ------------------------
# DESeq2 Analysis
# ------------------------
pdata_matched$mgs_level <- relevel(pdata_matched$mgs_level, ref = "1")

dds <- DESeqDataSetFromMatrix(countData = counts_matched,
                              colData = pdata_matched,
                              design = ~ mgs_level)

dds <- DESeq(dds)
res_mgs <- results(dds, contrast = c("mgs_level", "4", "1"))
norm_counts <- counts(dds, normalized = TRUE)



# ------------------------
# DEG Summary
# ------------------------
lfc_thresholds <- c(log2(1.5), log2(2))
fdr_threshold <- 0.05

filter_deg_summary <- function(lfc_cutoff, fdr_cutoff, res) {
  filtered <- res[!is.na(res$padj) & res$padj < fdr_cutoff & abs(res$log2FoldChange) >= lfc_cutoff, ]
  data.frame(
    LFC_Cutoff = lfc_cutoff,
    FDR_Cutoff = fdr_cutoff,
    Upregulated = sum(filtered$log2FoldChange >= lfc_cutoff),
    Downregulated = sum(filtered$log2FoldChange <= -lfc_cutoff),
    Total_DEGs = nrow(filtered)
  )
}

results_summary <- do.call(rbind, lapply(lfc_thresholds, filter_deg_summary, fdr_cutoff = fdr_threshold, res = res_mgs))
write.csv(results_summary, "results_all/tables/deg_summary_mgs_1_vs_4_all.csv")

# DEG Barplot
df_long <- pivot_longer(results_summary, cols = c("Upregulated", "Downregulated", "Total_DEGs"),
                        names_to = "Category", values_to = "Count")
df_long$LFC_Cutoff <- factor(df_long$LFC_Cutoff, labels = c("1.5 Fold", "2 Fold"))

png("results_all/plots/deg_counts_plot_mgs_1_vs_4_all.png", width = 800, height = 600)
ggplot(df_long, aes(x = LFC_Cutoff, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "DEG Counts: mgs_level 4 vs 1 (All Samples)",
       x = "Fold Change Threshold", y = "Number of DEGs") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal()
dev.off()

# ------------------------
# Save Results
# ------------------------
write.csv(as.data.frame(res_mgs), "results_all/tables/deseq2_results_mgs_4_vs_1_all.csv")

res_sig <- res_mgs[!is.na(res_mgs$padj) & res_mgs$padj < fdr_threshold & abs(res_mgs$log2FoldChange) > log2(1.5), ]
write.csv(as.data.frame(res_sig), "results_all/tables/deseq2_sig_results_mgs_1_vs_4_all.csv")
write.csv(norm_counts, "results_all/tables/normalized_counts_mgs_1_vs_4_all.csv")

# ------------------------
# Plots
# ------------------------
png("results_all/plots/ma_plot_mgs_1_vs_4_all.png", width = 800, height = 600)
plotMA(res_mgs, ylim = c(-4, 4), main = "MA Plot: mgs_level 4 vs 1 (All Samples)")
dev.off()

res_df <- as.data.frame(res_mgs)
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > log2(1.5), "Yes", "No")

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  labs(title = "Volcano Plot: mgs_level 4 vs 1 (All Samples)",
       x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme_minimal()

ggsave("results_all/plots/volcano_plot_mgs_1_vs_4_all.png", volcano_plot, width = 8, height = 6)
saveRDS(volcano_plot, file = "results_all/rds_plots/volcano_plot_mgs_1_vs_4_all.rds")

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
    labs(title = "Top 10 DEGs: mgs_level 4 vs 1 (All Samples)",
         x = "Gene", y = "Log2 Fold Change") +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
    theme_minimal()
  ggsave("results_all/plots/top_degs_plot_mgs_1_vs_4_all.png", top_degs_plot, width = 8, height = 6)
  saveRDS(top_degs_plot, file = "results_all/rds_plots/top_degs_plot_mgs_1_vs_4_all.rds")
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

# GO
go_dot <- dotplot(ego, showCategory = 20) + ggtitle("GO: Biological Process Enrichment")
png("results_all/plots/GO_BP_enrichment_1_vs_4.png", width = 1000, height = 800)
print(go_dot)
dev.off()
saveRDS(go_dot, file = "results_all/rds_plots/go_bp_dotplot_mgs_1_vs_4_all.rds")

# KEGG
kegg_dot <- dotplot(ekegg, showCategory = 20) + ggtitle("KEGG Pathway Enrichment")
png("results_all/plots/KEGG_enrichment_1_vs_4.png", width = 1000, height = 800)
print(kegg_dot)
dev.off()
saveRDS(kegg_dot, file = "results_all/rds_plots/kegg_dotplot_mgs_1_vs_4_all.rds")

# Reactome
react_bar <- barplot(react, showCategory = 20) + ggtitle("Reactome Pathway Enrichment")
png("results_all/plots/Reactome_enrichment_1_vs_4.png", width = 1000, height = 800)
print(react_bar)
dev.off()
saveRDS(react_bar, file = "results_all/rds_plots/reactome_barplot_mgs_1_vs_4_all.rds")

write.csv(as.data.frame(ego), "results_all/tables/go_bp_enrichment_1_vs_4.csv")
write.csv(as.data.frame(ekegg), "results_all/tables/kegg_enrichment_1_vs_4.csv")
write.csv(as.data.frame(react), "results_all/tables/reactome_enrichment_1_vs_4.csv")

# ------------------------
# Heatmap
# ------------------------
top_genes_heatmap <- if(nrow(res_sig) >= 30) rownames(head(res_sig[order(res_sig$padj), ], 30)) else rownames(res_sig)

if(length(top_genes_heatmap) > 0) {
  top_counts <- norm_counts[top_genes_heatmap, ]
  top_counts_scaled <- t(scale(t(top_counts)))
  
  annotation_col <- as.data.frame(colData(dds)[, "mgs_level", drop = FALSE])
  sample_order <- rownames(annotation_col)[order(annotation_col$mgs_level)]
  top_counts_sorted <- top_counts_scaled[, sample_order]
  annotation_col_sorted <- annotation_col[sample_order, , drop = FALSE]
  
  breaks <- seq(-2, 2, length.out = 51)
  
  pheatmap(top_counts_sorted,
           breaks = breaks,
           cluster_rows = TRUE, cluster_cols = FALSE,
           annotation_col = annotation_col_sorted,
           show_rownames = TRUE, show_colnames = TRUE,
           fontsize_row = 8, fontsize_col = 9,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
           main = paste("Heatmap of Top", length(top_genes_heatmap), "DEGs: mgs_level 4 vs 1 (All Samples)"),
           filename = "results_all/plots/heatmap_top_DEGs_mgs_1_vs_4_all.png")
  
  saveRDS(top_counts_sorted, "results_all/rds_plots/heatmap_matrix_top_DEGs_mgs_1_vs_4_all.rds")
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
  labs(title = "UMAP of Normalized Counts: mgs_level Groups (All Samples)",
       x = "UMAP1", y = "UMAP2") +
  theme_minimal() +
  scale_color_manual(values = c("1" = "blue", "4" = "red"))

ggsave("results_all/plots/umap_plot_mgs_1_vs_4_all.png", umap_plot, width = 8, height = 6)
saveRDS(umap_plot, file = "results_all/rds_plots/umap_plot_mgs_1_vs_4_all.rds")

# ------------------------
# Session Info
# ------------------------
writeLines(capture.output(sessionInfo()), "results_all/session_info_mgs_1_vs_4_all.txt")

cat("\nAnalysis complete for ALL samples!\n")
cat("Total samples analyzed:", nrow(pdata_matched), "\n")
cat("mgs_level group sizes:\n")
print(table(pdata_matched$mgs_level))
cat("\nSummary of significant DESeq2 results (res_sig):\n")
print(summary(res_sig))
