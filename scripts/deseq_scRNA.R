#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

count_file <- args[1]
meta_file  <- args[2]

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
})

# -----------------------------
# Read featureCounts output properly
# -----------------------------
counts_raw <- read.table(count_file,
                         header=TRUE,
                         row.names=1,
                         comment.char="#")

# Remove annotation columns
counts <- counts_raw[, 6:ncol(counts_raw)]

# Clean column names
colnames(counts) <- sub("_trimmed.bam$", "", colnames(counts))

# Remove ENSG version numbers
rownames(counts) <- sub("\\..*", "", rownames(counts))

meta <- read.table(meta_file, header=TRUE)
meta$sample <- as.character(meta$sample)

counts <- counts[, meta$sample]
rownames(meta) <- meta$sample

# Set Normal as reference
meta$condition <- factor(meta$condition,
                         levels = c("Normal","Tumor"))

# -----------------------------
# Convert ENSG → SYMBOL
# -----------------------------
gene_map <- bitr(rownames(counts),
                 fromType="ENSEMBL",
                 toType="SYMBOL",
                 OrgDb=org.Hs.eg.db)

gene_map <- gene_map[!is.na(gene_map$SYMBOL), ]

counts_df <- as.data.frame(counts)
counts_df$ENSEMBL <- rownames(counts_df)

counts_merged <- merge(counts_df, gene_map, by="ENSEMBL")

counts_final <- counts_merged %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), sum))

counts_final <- as.data.frame(counts_final)
rownames(counts_final) <- counts_final$SYMBOL
counts_final$SYMBOL <- NULL

counts <- as.matrix(counts_final)

# -----------------------------
# DESeq2
# -----------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta,
  design = ~ condition
)

dds <- DESeq(dds)
res <- results(dds)

# -----------------------------
# MA Plot
# -----------------------------
dir.create("plots", showWarnings = FALSE)

png("MA_plot.png", width=800, height=600)
plotMA(res, ylim=c(-5,5))
dev.off()

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

res_df$pvalue[is.na(res_df$pvalue)] <- 1
res_df$padj[is.na(res_df$padj)] <- 1

write.csv(res_df, "deseq_results.csv")

# -----------------------------
# PCA Plot
# -----------------------------
vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png("PCA_plot.png", width=900, height=700)

ggplot(pcaData, aes(PC1, PC2, color=condition, label=name)) +
  geom_point(size=4) +
  geom_text(vjust=1.5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot (Tumor vs Normal)") +
  theme_minimal()

dev.off()

# -----------------------------
# Volcano Plot
# -----------------------------
res_df$Significance <- "Not Significant"
res_df$Significance[res_df$pvalue < 0.05 & res_df$log2FoldChange > 1]  <- "Up"
res_df$Significance[res_df$pvalue < 0.05 & res_df$log2FoldChange < -1] <- "Down"

png("Volcano_plot.png", width=900, height=700)

ggplot(res_df, aes(x=log2FoldChange, y=-log10(pvalue), color=Significance)) +
  geom_point(alpha=0.6, size=1.5) +
  scale_color_manual(values=c("Up"="red",
                              "Down"="blue",
                              "Not Significant"="grey")) +
  geom_vline(xintercept=c(-1,1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  theme_minimal() +
  theme(text=element_text(size=12)) +
  labs(title="Volcano Plot (Tumor vs Normal)",
       x="Log2 Fold Change",
       y="-Log10 P-value")

dev.off()

# -----------------------------
# Heatmap (Top 20)
# -----------------------------
topgenes <- head(order(res$pvalue),20)
mat <- assay(vsd)[topgenes, ]
mat <- mat - rowMeans(mat)

png("Heatmap_top20.png", width=900, height=800)

pheatmap(mat,
         annotation_col = meta["condition"],
         show_rownames = TRUE,
         main="Top 20 Differentially Expressed Genes")

dev.off()

# -----------------------------
# Upregulated Genes
# -----------------------------
up <- res_df[res_df$log2FoldChange > 1 & res_df$pvalue < 0.05, ]
up <- up[order(up$log2FoldChange, decreasing=TRUE), ]
up_top <- head(up, 20)

png("Upregulated_genes.png", width=1000, height=900)

ggplot(up_top, aes(x=reorder(gene, log2FoldChange), y=log2FoldChange)) +
  geom_bar(stat="identity", fill="red") +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.y = element_text(size=12, color="black")) +
  labs(title="Top 20 Upregulated Genes (Tumor vs Normal)",
       x="Gene",
       y="Log2 Fold Change")

dev.off()

# -----------------------------
# Downregulated Genes
# -----------------------------
down <- res_df[res_df$log2FoldChange < -1 & res_df$pvalue < 0.05, ]
down <- down[order(down$log2FoldChange), ]
down_top <- head(down, 20)

png("Downregulated_genes.png", width=1000, height=900)

ggplot(down_top, aes(x=reorder(gene, log2FoldChange), y=log2FoldChange)) +
  geom_bar(stat="identity", fill="blue") +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.y = element_text(size=12, color="black")) +
  labs(title="Top 20 Downregulated Genes (Tumor vs Normal)",
       x="Gene",
       y="Log2 Fold Change")

dev.off()
# -----------------------------
# Enrichment Analysis
# -----------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

sig <- res_df[res_df$pvalue < 0.05 & abs(res_df$log2FoldChange) > 1, ]

gene_symbols <- sig$gene

gene_ids <- bitr(gene_symbols,
                 fromType="SYMBOL",
                 toType="ENTREZID",
                 OrgDb=org.Hs.eg.db)

ego <- enrichGO(gene = gene_ids$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                readable = TRUE)

png("GO_enrichment.png", width=1000, height=900)
dotplot(ego, showCategory=20)
dev.off()
# -----------------------------
# KEGG Enrichment
# -----------------------------
ekegg <- enrichKEGG(
  gene = gene_ids$ENTREZID,
  organism = "hsa"
)

png("KEGG_enrichment.png", width=1000, height=900)

dotplot(ekegg, showCategory=20) +
  ggtitle("KEGG Pathway Enrichment")

dev.off()
# -----------------------------
# Save KEGG pathway links
# -----------------------------
kegg_df <- as.data.frame(ekegg)

dir.create("results/deseq", showWarnings = FALSE, recursive = TRUE)

if(nrow(kegg_df) > 0){

  kegg_links <- paste0(
    "https://www.kegg.jp/pathway/",
    kegg_df$ID
  )

  kegg_table <- data.frame(
    Pathway = kegg_df$Description,
    Link = kegg_links
  )

  write.csv(
    kegg_table,
    "results/deseq/kegg_links.csv",
    row.names = FALSE
  )

} else {

  write.csv(
    data.frame(Pathway="No enriched KEGG pathways",
               Link=""),
    "results/deseq/kegg_links.csv",
    row.names = FALSE
  )

}
# -----------------------------
# Move plots to results folder
# -----------------------------
plots_folder <- file.path("results", "plots")
dir.create(plots_folder, showWarnings = FALSE, recursive = TRUE)

png_plots <- c(
  "MA_plot.png",
  "PCA_plot.png",
  "Volcano_plot.png",
  "Heatmap_top20.png",
  "Upregulated_genes.png",
  "Downregulated_genes.png",
  "GO_enrichment.png",
  "KEGG_enrichment.png"
)
for(p in png_plots){
  if(file.exists(p)){
    file.rename(p, file.path(plots_folder, p))
  }
}

cat("All DESeq2 plots successfully moved to results/plots\n")
