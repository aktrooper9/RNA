library(sva)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tibble)
library(readr)
library(limma)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ReactomePA)
library(cowplot)

counts_data <- read.table("merged_all_counts.tsv", header = TRUE, row.names = 1, sep = "\t")
metadata <- read.table("combined_geo_accession_study_condition.tsv", header = TRUE, sep = "\t")
rownames(metadata) <- metadata$sample_geo_accession

metadata$condition <- as.factor(metadata$condition)
metadata$study <- as.factor(metadata$study)

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = metadata,
                              design = ~ condition)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)
res <- results(dds)
write.csv(as.data.frame(res), "DESeq2_results.csv")
cat("DESeq2 Results Summary:\n")
summary(res)

vsd <- vst(dds, blind = FALSE)
pca_plot <- plotPCA(vsd, intgroup = c("condition", "study")) + ggtitle("PCA Before Batch Correction")
ggsave("pca_before_batch_correction.png", pca_plot, width = 8, height = 6)



sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
up_genes <- subset(sig_genes, log2FoldChange > 1)
down_genes <- subset(sig_genes, log2FoldChange < -1)
write.csv(as.data.frame(sig_genes), "significant_genes.csv")
cat("Number of significantly differentially expressed genes:", nrow(sig_genes), "\n")
cat("Upregulated genes:", nrow(up_genes), "\n")
cat("Downregulated genes:", nrow(down_genes), "\n")
if (nrow(sig_genes) == 0) warning("No significant genes found with current thresholds.")




norm_counts <- counts(dds, normalized = TRUE)
gene_vector <- rownames(sig_genes)
  heatmap_data <- norm_counts[gene_vector, ]
  pheatmap(heatmap_data, 
           annotation_col = metadata[, c("condition", "study")],
           color = colorRampPalette(brewer.pal(9, "RdBu"))(100),
           main = "Heatmap of Significant Genes",
           filename = "heatmap_significant_genes.png",
           width = 10, height = 8)
  
  
  up_entrez <- if (nrow(up_genes) > 0) rownames(up_genes) else character(0)
  down_entrez <- if (nrow(down_genes) > 0) rownames(down_genes) else character(0)
  if (length(up_entrez) == 0) cat("No upregulated genes for enrichment analysis.\n")
  if (length(down_entrez) == 0) cat("No downregulated genes for enrichment analysis.\n")
  
  
  
  
  
  
  
  
    ego_up <- enrichGO(gene = up_entrez,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2,
                       readable = TRUE)
    ego_up <- simplify(ego_up, cutoff = 0.7, by = "p.adjust", select_fun = min)
    write.csv(as.data.frame(ego_up), "go_enrichment_up.csv")
    dotplot_up <- dotplot(ego_up, showCategory = 15) + ggtitle("Top 15 GO Biological Processes (Upregulated)")
    ggsave("go_dotplot_up.png", dotplot_up, width = 10, height = 8)
    barplot_up <- barplot(ego_up, showCategory = 15) + ggtitle("Top 15 GO Biological Processes (Upregulated)")
    ggsave("go_barplot_up.png", barplot_up, width = 10, height = 8)
  
  
    ego_down <- enrichGO(gene = down_entrez,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID",
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2,
                         readable = TRUE)
    ego_down <- simplify(ego_down, cutoff = 0.7, by = "p.adjust", select_fun = min)
    write.csv(as.data.frame(ego_down), "go_enrichment_down.csv")
    dotplot_down <- dotplot(ego_down, showCategory = 15) + ggtitle("Top 15 GO Biological Processes (Downregulated)")
    ggsave("go_dotplot_down.png", dotplot_down, width = 10, height = 8)
    barplot_up <- barplot(ego_down, showCategory = 15) + ggtitle("Top 15 GO Biological Processes (Downregulated)")
    ggsave("go_barplot_down.png", barplot_up, width = 10, height = 8)
    
    
    
    
    
    
      kegg_up <- enrichKEGG(gene = up_entrez, organism = "hsa", pvalueCutoff = 0.05)
      write.csv(as.data.frame(kegg_up), "kegg_enrichment_up.csv")
      kegg_dotplot_up <- dotplot(kegg_up, showCategory = 15) + ggtitle("KEGG Pathways (Upregulated)")
      ggsave("kegg_dotplot_up.png", kegg_dotplot_up, width = 10, height = 8)
    
    
      kegg_down <- enrichKEGG(gene = down_entrez, organism = "hsa", pvalueCutoff = 0.05)
      write.csv(as.data.frame(kegg_down), "kegg_enrichment_down.csv")
      kegg_dotplot_down <- dotplot(kegg_down, showCategory = 15) + ggtitle("KEGG Pathways (Downregulated)")
      ggsave("kegg_dotplot_down.png", kegg_dotplot_down, width = 10, height = 8)
    
    
    # Reactome enrichment analysis
      reactome_res <- enrichPathway(gene = up_entrez, pvalueCutoff = 0.05, readable = TRUE)
      write.csv(as.data.frame(reactome_res), "reactome_enrichment.csv")
      reactome_dotplot <- dotplot(reactome_res, showCategory = 15) + ggtitle("Reactome Pathway Enrichment")
      ggsave("reactome_dotplot.png", reactome_dotplot, width = 10, height = 8)
      cnetplot_kegg <- cnetplot(kegg_up, showCategory = 5) + ggtitle("KEGG Pathway Network")
      ggsave("kegg_cnetplot.png", cnetplot_kegg, width = 10, height = 8)
    
    
    # Reactome enrichment analysis
      reactome_res <- enrichPathway(gene = down_entrez, pvalueCutoff = 0.05, readable = TRUE)
      write.csv(as.data.frame(reactome_res), "reactome_enrichment.csv")
      reactome_dotplot <- dotplot(reactome_res, showCategory = 15) + ggtitle("Reactome Pathway Enrichment")
      ggsave("reactome_down_dotplot.png", reactome_dotplot, width = 10, height = 8)
      cnetplot_kegg <- cnetplot(kegg_up, showCategory = 5) + ggtitle("KEGG Pathway Network")
      ggsave("kegg_down_cnetplot.png", cnetplot_kegg, width = 10, height = 8)
    
    
    # GSEA analysis
    ranked_gene_vector <- res$log2FoldChange
    names(ranked_gene_vector) <- rownames(res)
    ranked_gene_vector <- sort(ranked_gene_vector, decreasing = TRUE)
    gsea_res <- gseGO(geneList = ranked_gene_vector, 
                      OrgDb = org.Hs.eg.db,
                      ont = "BP", 
                      pvalueCutoff = 0.05)
    write.csv(as.data.frame(gsea_res), "gsea_results.csv")
    gsea_dotplot <- dotplot(gsea_res, showCategory = 15) + ggtitle("GSEA: Top 15 GO Biological Processes")
    ggsave("gsea_dotplot.png", gsea_dotplot, width = 10, height = 8)
  
  
