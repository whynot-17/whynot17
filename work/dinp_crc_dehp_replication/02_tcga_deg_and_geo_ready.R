#!/usr/bin/env Rscript

# Run after 01_preprocess_crc_transcriptomes.R
# Reproduces the paper's TCGA CRC-vs-normal limma screen:
# nominal P < 0.05 and |log2FC| > 1.
# GEO merged data are exported as the machine-learning/validation matrix.

options(stringsAsFactors = FALSE)
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
outdir <- file.path(root, "outputs", "dinp_crc_dehp_replication")

if (!requireNamespace("limma", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install("limma", ask = FALSE, update = FALSE)
}
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table", repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({library(limma); library(data.table)})

read_expr <- function(path) {
  x <- fread(path, data.table = FALSE, check.names = FALSE)
  rn <- x[[1]]; x[[1]] <- NULL
  m <- as.matrix(x); storage.mode(m) <- "numeric"; rownames(m) <- rn
  m
}

tcga <- read_expr(file.path(outdir, "tcga_coread_expression_gene.tsv.gz"))
ph <- fread(file.path(outdir, "tcga_coread_pheno.tsv"), data.table = FALSE)
rownames(ph) <- ph$sample_id
ph <- ph[colnames(tcga), , drop = FALSE]
ph$group <- factor(ph$group, levels = c("Normal", "Tumor"))

design <- model.matrix(~ 0 + group, data = ph)
colnames(design) <- levels(ph$group)
fit <- lmFit(tcga, design)
cont <- makeContrasts(Tumor - Normal, levels = design)
fit2 <- eBayes(contrasts.fit(fit, cont))
res <- topTable(fit2, number = Inf, sort.by = "P")
res$gene <- rownames(res)
res$paper_hit <- res$P.Value < 0.05 & abs(res$logFC) > 1
res$fdr_hit <- res$adj.P.Val < 0.05 & abs(res$logFC) > 1
res <- res[, c("gene", setdiff(names(res), "gene"))]
fwrite(res, file.path(outdir, "TCGA_COADREAD_limma_all.tsv"), sep = "\t", quote = FALSE)
fwrite(res[res$paper_hit, ], file.path(outdir, "TCGA_COADREAD_DEG_paper_threshold.tsv"), sep = "\t", quote = FALSE)
fwrite(res[res$fdr_hit, ], file.path(outdir, "TCGA_COADREAD_DEG_FDR_threshold.tsv"), sep = "\t", quote = FALSE)

# Machine-learning/validation matrix: use ComBat version by default, retain uncorrected from step 01.
geo <- read_expr(file.path(outdir, "GEO_merged_combat.tsv.gz"))
gph <- fread(file.path(outdir, "GEO_merged_pheno.tsv"), data.table = FALSE)
rownames(gph) <- gph$sample_id
gph <- gph[colnames(geo), , drop = FALSE]

# samples x genes, ready for LASSO/SVM/RF after target intersection/PPI candidate selection.
ml <- data.frame(sample_id = colnames(geo), group = gph$group, t(geo), check.names = FALSE)
fwrite(ml, file.path(outdir, "GEO_ML_ready_samples_x_genes.tsv.gz"), sep = "\t", quote = FALSE)

cat("TCGA limma done. Paper-threshold DEGs:", sum(res$paper_hit), "\n")
cat("FDR-threshold DEGs:", sum(res$fdr_hit), "\n")
cat("GEO ML matrix:", nrow(ml), "samples x", ncol(ml)-2, "genes\n")
