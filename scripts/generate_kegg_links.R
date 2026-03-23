#!/usr/bin/env Rscript

library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(pathview)

# --------------------------------------------------
# Locate DESeq2 results inside Nextflow work folders
# --------------------------------------------------

cat("Searching for DESeq2 results inside Nextflow work directory...\n")

files <- list.files(
  path = "../../",
  pattern = "deseq_results.csv",
  recursive = TRUE,
  full.names = TRUE
)

if(length(files) == 0){
  stop("ERROR: deseq_results.csv not found in pipeline work folders")
}

res_file <- files[1]

cat("DESeq2 results found at:", res_file, "\n")

res <- read.csv(res_file, row.names = 1)

# --------------------------------------------------
# Check required columns
# --------------------------------------------------

required_cols <- c("log2FoldChange","padj")

if(!all(required_cols %in% colnames(res))){
  stop("ERROR: Required DESeq2 columns missing")
}

# --------------------------------------------------
# Filter significant genes
# --------------------------------------------------

sig_genes <- rownames(res)[res$padj < 0.05]
sig_genes <- sig_genes[!is.na(sig_genes)]

cat("Significant genes found:", length(sig_genes), "\n")

if(length(sig_genes) == 0){
  stop("No significant genes found")
}

# --------------------------------------------------
# Convert SYMBOL → ENTREZ
# --------------------------------------------------

gene_df <- bitr(sig_genes,
                fromType="SYMBOL",
                toType="ENTREZID",
                OrgDb=org.Hs.eg.db)

if(nrow(gene_df) == 0){
  stop("Gene conversion failed")
}

cat("Gene ID conversion successful\n")

# --------------------------------------------------
# KEGG enrichment
# --------------------------------------------------

kegg_res <- enrichKEGG(
  gene = gene_df$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

if(is.null(kegg_res)){
  stop("KEGG enrichment failed")
}

kegg_df <- as.data.frame(kegg_res)

if(nrow(kegg_df) == 0){
  stop("No KEGG pathways enriched")
}

cat("KEGG pathways identified:", nrow(kegg_df), "\n")

# --------------------------------------------------
# Prepare fold change vector for Pathview
# --------------------------------------------------

gene_fc <- res$log2FoldChange
names(gene_fc) <- rownames(res)

gene_map <- bitr(names(gene_fc),
                 fromType="SYMBOL",
                 toType="ENTREZID",
                 OrgDb=org.Hs.eg.db)

gene_fc <- gene_fc[gene_map$SYMBOL]
names(gene_fc) <- gene_map$ENTREZID

# --------------------------------------------------
# Create results folders
# --------------------------------------------------

dir.create("results/deseq", showWarnings=FALSE, recursive=TRUE)
dir.create("results/pathview", showWarnings=FALSE, recursive=TRUE)

# --------------------------------------------------
# Generate pathway diagrams
# --------------------------------------------------

cat("Generating pathway diagrams using Pathview...\n")

top_paths <- head(kegg_df$ID, 5)

for(pid in top_paths){

  tryCatch({

    pathview(
      gene.data = gene_fc,
      pathway.id = pid,
      species = "hsa",
      out.suffix = pid,
      kegg.native = TRUE
    )

  }, error=function(e){
    cat("Pathview failed for:", pid, "\n")
  })

}

# --------------------------------------------------
# Create clickable KEGG links
# --------------------------------------------------

kegg_links <- data.frame(
  Pathway = kegg_df$Description,
  KEGG_ID = kegg_df$ID,
  Link = paste0("https://www.kegg.jp/pathway/", kegg_df$ID)
)

# Save CSV required by Nextflow
write.csv(
  kegg_links,
  "kegg_links.csv",
  row.names = FALSE
)

# Optional copy to results folder
dir.create("results/deseq", showWarnings = FALSE, recursive = TRUE)

write.csv(
  kegg_links,
  "results/deseq/kegg_links.csv",
  row.names = FALSE
)

cat("KEGG links saved\n")
cat("KEGG analysis completed successfully\n")
