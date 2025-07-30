library(DESeq2)
library(pheatmap)
library(ggplot2)
library(edgeR)
#library(clusterProfiler)
#library(org.Hs.eg.db)
#library(enrichplot)


# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE115828", "file=GSE115828_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")

# load gene annotations from GEO
apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
rownames(annot) <- annot$GeneID

# load counts from local
#path <- ("counts.tsv")
#tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")

# load gene annotations from local
#apath <- ("annot.tsv")
#annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
#rownames(annot) <- annot$GeneID

# sample selection
gsms <- paste0("01011102223501011424330311122410021112114133134040",
               "02412203111210350340011014513413042141110101004111",
               "11001050301102110005000113104300103221202301051510",
               "44011041041120040230001310013440011512044010005011",
               "21303104432315103021500433101403434444433444011005",
               "10012311101000011444433443333100010012110000101113",
               "3314100154040231100X334511300102344000302133325024",
               "21101311203000010150110011000411015311500125401011",
               "31420012110511043105011111330403150331040400400511",
               "11012243001202114100110111103501101011100124010001",
               "041011000430110220014")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
tbl <- tbl[ ,sel]

# group membership for samples
gs <- factor(sml)
groups <- make.names(c("M23","F23","F4","F1","M1","M4"))
levels(gs) <- groups
sample_info <- data.frame(Group = gs, row.names = colnames(tbl))

# pre-filter low count genes
# keep genes with at least N counts > 10, where N = size of smallest group
keep <- rowSums( tbl >= 10 ) >= min(table(gs))
tbl <- tbl[keep, ]

ds <- DESeqDataSetFromMatrix(countData=tbl, colData=sample_info, design= ~Group)

cts <- list(c("Group",groups[1],groups[2]),
            c("Group",groups[3],groups[4]),
            c("Group",groups[3],groups[6]),
            c("Group",groups[4],groups[5]),
            c("Group",groups[5],groups[6]))  # contrasts of interest

ds <- DESeq(ds, test="LRT", reduced = ~ 1)  # Use LRT for all-around gene ranking

# extract results for top genes table
r <- results (ds, alpha=0.05, pAdjustMethod ="fdr")
summary(r)
r_filtered <- r[abs(r$log2FoldChange)>1.5 & r$padj <0.05, ]
summary(r_filtered)

tT <- r[order(r$padj)[1:250],]
tT <- merge(as.data.frame(tT), annot, by=0, sort=F)

tT <- subset(tT, select=c("GeneID","padj","pvalue","stat","baseMean","Symbol","Description"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

plotDispEsts(ds, main="GSE115828 Dispersion Estimates")

# create histogram plot of p-values
hist(r$padj, breaks=seq(0, 1, length = 21), col = "grey", border = "white",
     xlab = "", ylab = "", main = "GSE115828 Frequencies of padj-values")

# Wald test to obtain contrast-specific results
#ds <- DESeq(ds, test="Wald", sfType="poscount")
#r <- results (ds, contrast=cts[[1]], alpha=0.05, pAdjustMethod = "fdr")

vsd <- vst(ds, blind =TRUE)

for (i in seq_along(cts)) {
  # Define contrast name early for messaging and filenames
  contrast_name <- paste(cts[[i]][2], "vs", cts[[i]][3])
  
  # Perform DESeq2 analysis for each contrast
  r <- results(ds, contrast = cts[[i]])
  
  # Filter results: adjusted p-value < 0.05 and |LFC| > 1.5
  r_filtered <- r[r$padj < 0.05 & abs(r$log2FoldChange) > 1.5, ]
  
  # Check if there are significant genes
  if (nrow(r_filtered) == 0) {
    cat("No significant genes for contrast:", contrast_name, "\n")
    next
  }
  
  # Find common genes between filtered results and vst assay rows
  common_genes <- intersect(rownames(r_filtered), rownames(assay(vsd)))
  if (length(common_genes) == 0) {
    cat("No overlapping genes for contrast:", contrast_name, "\n")
    next
  }
  
  # Convert filtered results to data frame and subset to common genes
  r_filtered_df <- as.data.frame(r_filtered)[common_genes, , drop = FALSE]
  
  # Extract gene annotation info for these genes
  gene_info <- rowData(ds)[common_genes, ]
  
  # Combine filtered results with gene info
  result_with_info <- cbind(r_filtered_df, gene_info)
  
  # Save significant genes to CSV
  significant_genes_file <- paste0("significant_genes_", contrast_name, ".csv")
  write.csv(result_with_info, file = significant_genes_file)
  
  # Get VST-transformed counts for common genes
  vst_counts <- assay(vsd)[common_genes, ]
  
  # Select samples belonging to the two groups in the contrast
  groups_to_include <- c(cts[[i]][2], cts[[i]][3])
  coldata <- colData(ds)
  selected_samples <- coldata[coldata$Group %in% groups_to_include, ]
  
  # Sort samples by group (group1 then group2)
  selected_samples <- selected_samples[order(selected_samples$Group), ]
  
  # Subset vst_counts columns to these selected samples in sorted order
  vst_counts <- vst_counts[, rownames(selected_samples)]
  
  # Annotation for heatmap columns
  annotation_col <- data.frame(Group = selected_samples$Group)
  rownames(annotation_col) <- colnames(vst_counts)
  
  # Replace rownames of vst_counts with gene symbols for better heatmap labels
  if (any(is.na(result_with_info$Symbol))) {
    warning("NA values found in gene symbols; replacing with gene IDs.")
    result_with_info$Symbol[is.na(result_with_info$Symbol)] <- rownames(result_with_info)[is.na(result_with_info$Symbol)]
  }
  rownames(vst_counts) <- result_with_info$Symbol
  
  # Save heatmap to PNG file
  png_filename <- paste0("heatmap_", contrast_name, ".png")
  png(png_filename, width = 1000, height = 800)
  pheatmap(vst_counts,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           scale = "row",
           main = paste("Heatmap for", contrast_name),
           show_rownames = TRUE,
           show_colnames = TRUE,
           annotation_col = annotation_col,
           fontsize = 12)
  dev.off()
  
  # Volcano plot using full results 'r' (not filtered)
  volcano_plot <- ggplot(as.data.frame(r), aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1.5), alpha = 0.5) +
    scale_color_manual(values = c("black", "red")) +
    theme_minimal() +
    labs(title = paste("Volcano Plot for", contrast_name),
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-Value") +
    theme(legend.position = "none")
  ggsave(paste0("volcano_", contrast_name, ".png"), plot = volcano_plot, width = 8, height = 6)
  
  # MD plot: mean normalized counts vs log2FoldChange
  r$mean_count <- rowMeans(counts(ds, normalized = TRUE)[rownames(r), ])
  md_plot <- ggplot(as.data.frame(r), aes(x = mean_count, y = log2FoldChange)) +
    geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1.5), alpha = 0.5) +
    scale_color_manual(values = c("black", "red")) +
    theme_minimal() +
    labs(title = paste("MD Plot for", contrast_name),
         x = "Mean Normalized Count",
         y = "Log2 Fold Change") +
    theme(legend.position = "none")
  ggsave(paste0("md_plot_", contrast_name, ".png"), plot = md_plot, width = 8, height = 6)
  
  # PCA plot using filtered vst counts
  pca_data <- prcomp(t(vst_counts))
  pca_df <- data.frame(pca_data$x, Group = selected_samples$Group)
  pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 3) +
    labs(title = paste("PCA Plot for", contrast_name),
         x = "PC1",
         y = "PC2") +
    theme_minimal() +
    scale_color_manual(values = c('red','blue','green','purple','orange'))
  ggsave(paste0("pca_", contrast_name, ".png"), plot = pca_plot, width = 8, height = 6)
  
  # Inform about saved files
  cat("Saved heatmap, volcano plot, MD plot, PCA plot, and significant genes file for", contrast_name, "\n")
}
``
