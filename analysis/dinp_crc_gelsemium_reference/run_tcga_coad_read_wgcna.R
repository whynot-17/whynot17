#!/usr/bin/env Rscript

# Full-transcriptome WGCNA for TCGA COAD + READ primary tumours.
# The fresh 41-gene DINP–CRC set is overlaid only after network construction.

suppressPackageStartupMessages({
  library(WGCNA)
  library(dynamicTreeCut)
  library(fastcluster)
})
options(stringsAsFactors = FALSE)
# WGCNA in the installed R build requires at least two threads. Keep the
# execution deliberately small and reproducible on the local workstation.
allowWGCNAThreads(nThreads = 2)

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")) else getwd()
ROOT <- normalizePath(file.path(script_dir, "../.."), winslash = "/", mustWork = TRUE)
RUNTIME <- Sys.getenv("WGCNA_RUNTIME", unset = "D:/whynot17/work/tcga_coad_read_wgcna_runtime")
OUTPUT_SUBDIR <- Sys.getenv("WGCNA_OUTPUT_SUBDIR", unset = "tcga_coad_read_wgcna")
NETWORK_TYPE <- Sys.getenv("WGCNA_NETWORK_TYPE", unset = "signed")
COR_TYPE <- Sys.getenv("WGCNA_COR_TYPE", unset = "bicor")
OUTPUT <- file.path(ROOT, "analysis", "dinp_crc_gelsemium_reference", "outputs", OUTPUT_SUBDIR)
INPUT_MATRIX <- file.path(RUNTIME, "tcga_coad_read_expression_top_variable.csv.gz")
INPUT_SELECTION <- file.path(OUTPUT, "tcga_coad_read_wgcna_selected_genes.csv")
INPUT_META <- file.path(OUTPUT, "tcga_coad_read_wgcna_sample_manifest.csv")
TARGET_FILE <- file.path(ROOT, "analysis", "dinp_crc_gelsemium_reference", "outputs", "query_41_genes.csv")
COAD_EXPR <- file.path(ROOT, "analysis", "dinp_crc_gelsemium_reference", "outputs", "tcga_coad_bulk_41_genes", "tcga_coad_expression_41_genes.csv")
READ_EXPR <- file.path(ROOT, "analysis", "dinp_crc_gelsemium_reference", "outputs", "tcga_read_bulk_41_genes", "tcga_read_expression_41_genes.csv")
dir.create(OUTPUT, recursive = TRUE, showWarnings = FALSE)

required <- c(INPUT_MATRIX, INPUT_SELECTION, INPUT_META, TARGET_FILE, COAD_EXPR, READ_EXPR)
if (any(!file.exists(required))) stop("Missing required input: ", paste(required[!file.exists(required)], collapse = "; "))
write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, quote = TRUE, na = "")
bh <- function(p) p.adjust(as.numeric(p), method = "BH")
safe_cor <- function(x, y) {
  if (sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0) return(NA_real_)
  cor(x, y, use = "pairwise.complete.obs", method = "pearson")
}

message("Reading transcriptome input")
raw <- read.csv(gzfile(INPUT_MATRIX), check.names = FALSE, stringsAsFactors = FALSE)
sample_ids <- raw[[1]]
datExpr <- as.data.frame(lapply(raw[-1], as.numeric), check.names = FALSE)
rownames(datExpr) <- sample_ids
if (anyNA(datExpr)) stop("The top-variable expression matrix contains NA values")
if (anyDuplicated(colnames(datExpr))) stop("The WGCNA input contains duplicate gene symbols")
metadata <- read.csv(INPUT_META, stringsAsFactors = FALSE)
metadata <- metadata[match(rownames(datExpr), metadata$sample_id), , drop = FALSE]
if (anyNA(metadata$sample_id) || !identical(rownames(datExpr), metadata$sample_id)) stop("Sample metadata does not match expression rows")
if (any(table(metadata$cohort) < 3)) stop("Each cohort needs at least three samples")

gsg <- goodSamplesGenes(datExpr, verbose = 2)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  metadata <- metadata[match(rownames(datExpr), metadata$sample_id), , drop = FALSE]
}
sample_qc <- data.frame(
  sample_id = rownames(datExpr), cohort = metadata$cohort,
  mean_expression = rowMeans(datExpr), sd_expression = apply(datExpr, 1, sd),
  missing_values = rowSums(is.na(datExpr))
)
write_csv(sample_qc, file.path(OUTPUT, "tcga_coad_read_wgcna_sample_qc.csv"))
sample_tree <- hclust(dist(datExpr), method = "average")
pdf(file.path(OUTPUT, "tcga_coad_read_wgcna_sample_dendrogram.pdf"), width = 12, height = 7)
plot(sample_tree, main = "TCGA COAD + READ primary-tumour clustering", xlab = "", sub = "No automatic sample removal")
dev.off()

cor_options <- if (COR_TYPE == "bicor") list(maxPOutliers = 0.05) else list(use = "p")
cor_fnc <- if (COR_TYPE == "pearson") "cor" else COR_TYPE
message("Selecting ", NETWORK_TYPE, "-network soft-threshold power with ", COR_TYPE)
powers <- c(seq(1, 10, by = 1), seq(12, 30, by = 2))
sft <- suppressWarnings(pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = NETWORK_TYPE,
  corFnc = cor_fnc,
  corOptions = cor_options,
  verbose = 2
))
sft_table <- sft$fitIndices
write_csv(sft_table, file.path(OUTPUT, "tcga_coad_read_wgcna_soft_threshold.csv"))
pdf(file.path(OUTPUT, "tcga_coad_read_wgcna_soft_threshold.pdf"), width = 10, height = 7)
plot(sft_table$Power, sft_table$SFT.R.sq, type = "b", pch = 19,
     xlab = "Soft-thresholding power", ylab = "Signed scale-free topology fit (R²)",
     main = "Soft-threshold selection")
abline(h = 0.80, col = "firebrick", lty = 2)
dev.off()
eligible_power <- sft_table$Power[is.finite(sft_table$SFT.R.sq) &
                                    sft_table$SFT.R.sq >= 0.80 &
                                    is.finite(sft_table$mean.k.) & sft_table$mean.k. > 1]
if (length(eligible_power)) {
  soft_power <- min(eligible_power)
  power_rule <- paste0("first power with ", NETWORK_TYPE, " R2 >= 0.80 and mean connectivity > 1")
} else {
  finite_fit <- which(is.finite(sft_table$SFT.R.sq))
  if (!length(finite_fit)) stop("No finite soft-threshold fit was returned")
  best_fit <- finite_fit[which.max(sft_table$SFT.R.sq[finite_fit])]
  soft_power <- sft_table$Power[best_fit]
  power_rule <- paste0("fallback: power with maximum finite ", NETWORK_TYPE, " R2")
}
message("Selected soft-threshold power: ", soft_power, " (", power_rule, ")")

message("Building ", NETWORK_TYPE, " WGCNA modules")
net <- blockwiseModules(
  datExpr, power = soft_power, TOMType = NETWORK_TYPE, networkType = NETWORK_TYPE,
  corType = COR_TYPE, maxPOutliers = 0.05, minModuleSize = 30,
  mergeCutHeight = 0.25, deepSplit = 2, pamRespectsDendro = FALSE,
  numericLabels = TRUE, saveTOMs = FALSE, maxBlockSize = ncol(datExpr), verbose = 2
)
module_colors <- labels2colors(net$colors)
names(module_colors) <- colnames(datExpr)
MEs <- orderMEs(net$MEs)
module_labels <- as.integer(sub("^ME", "", colnames(MEs)))
module_names <- labels2colors(module_labels)

pdf(file.path(OUTPUT, "tcga_coad_read_wgcna_module_dendrogram.pdf"), width = 14, height = 8)
for (block in seq_along(net$dendrograms)) {
  plot(net$dendrograms[[block]], main = paste("WGCNA gene dendrogram, block", block), xlab = "", sub = "Dynamic tree cut; merged signed modules")
  block_colors <- module_colors[net$blockGenes[[block]]]
  if (length(block_colors)) plotColorUnderTree(net$dendrograms[[block]], block_colors, main = "Module colors")
}
dev.off()

kme <- as.data.frame(cor(datExpr, MEs, use = "pairwise.complete.obs", method = "pearson"))
colnames(kme) <- paste0("kME_", module_names)
gene_membership <- data.frame(gene_symbol = colnames(datExpr), module = unname(module_colors[colnames(datExpr)]), network_input = TRUE)
gene_membership <- cbind(gene_membership, kme[gene_membership$gene_symbol, , drop = FALSE])
kme_cols <- grep("^kME_", colnames(gene_membership), value = TRUE)
gene_membership$max_abs_kME <- apply(abs(gene_membership[, kme_cols, drop = FALSE]), 1, max, na.rm = TRUE)
gene_membership$best_kME_module <- sub("^kME_", "", kme_cols[max.col(abs(gene_membership[, kme_cols, drop = FALSE]), ties.method = "first")])
write_csv(gene_membership, file.path(OUTPUT, "tcga_coad_read_wgcna_gene_module_membership.csv"))

target_genes <- unique(trimws(as.character(read.csv(TARGET_FILE, stringsAsFactors = FALSE)$gene_symbol)))
if (length(target_genes) != 41) stop("Expected exactly 41 target genes")
coad_target <- read.csv(COAD_EXPR, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)[, target_genes, drop = FALSE]
read_target <- read.csv(READ_EXPR, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)[, target_genes, drop = FALSE]
target_expr <- rbind(coad_target, read_target)
target_expr <- target_expr[match(rownames(datExpr), rownames(target_expr)), target_genes, drop = FALSE]
if (anyNA(target_expr)) stop("Target expression did not match all WGCNA samples")

target_z <- scale(target_expr)
target_score <- rowMeans(target_z, na.rm = TRUE)
cohort_numeric <- ifelse(metadata$cohort == "COAD", 1, 0)
module_trait <- do.call(rbind, lapply(seq_len(ncol(MEs)), function(i) {
  target_cor <- safe_cor(MEs[, i], target_score)
  cohort_cor <- safe_cor(MEs[, i], cohort_numeric)
  data.frame(module = module_names[i], module_size = sum(module_colors == module_names[i]),
             target_score_cor = target_cor, target_score_p = if (is.finite(target_cor)) corPvalueStudent(target_cor, nrow(datExpr)) else NA_real_,
             cohort_COAD_cor = cohort_cor, cohort_COAD_p = if (is.finite(cohort_cor)) corPvalueStudent(cohort_cor, nrow(datExpr)) else NA_real_)
}))
module_trait$target_score_BH_FDR <- bh(module_trait$target_score_p)
module_trait$cohort_COAD_BH_FDR <- bh(module_trait$cohort_COAD_p)
write_csv(module_trait, file.path(OUTPUT, "tcga_coad_read_wgcna_module_trait_correlations.csv"))

network_targets <- intersect(target_genes, colnames(datExpr))
target_overlay <- do.call(rbind, lapply(target_genes, function(gene) {
  in_network <- gene %in% colnames(datExpr)
  projected_cor <- vapply(seq_len(ncol(MEs)), function(i) safe_cor(as.numeric(target_expr[, gene]), MEs[, i]), numeric(1))
  best <- if (all(!is.finite(projected_cor))) NA_integer_ else which.max(abs(projected_cor))
  actual <- if (in_network) unname(module_colors[gene]) else NA_character_
  actual_kme_col <- if (isTRUE(in_network) && !is.na(actual) && actual != "grey" && paste0("kME_", actual) %in% colnames(kme)) paste0("kME_", actual) else NA_character_
  actual_kme <- if (!is.na(actual_kme_col)) unname(kme[gene, actual_kme_col]) else NA_real_
  data.frame(gene_symbol = gene, in_network_input = in_network, actual_module = actual,
             actual_kME = actual_kme,
             projected_module = if (is.na(best)) NA_character_ else module_names[best],
             projected_kME = if (is.na(best)) NA_real_ else projected_cor[best],
             max_abs_projected_kME = if (all(!is.finite(projected_cor))) NA_real_ else max(abs(projected_cor), na.rm = TRUE))
}))
write_csv(target_overlay, file.path(OUTPUT, "tcga_coad_read_wgcna_41_gene_overlay.csv"))

enrichment <- do.call(rbind, lapply(sort(unique(module_colors)), function(module) {
  in_module <- module_colors == module
  a <- sum(names(module_colors)[in_module] %in% network_targets)
  b <- sum(names(module_colors)[!in_module] %in% network_targets)
  c <- sum(in_module) - a
  d <- sum(!in_module) - b
  ft <- fisher.test(matrix(c(a, c, b, d), nrow = 2, byrow = TRUE), alternative = "greater")
  data.frame(module = module, module_size = sum(in_module), target_genes_in_network = length(network_targets), target_overlap = a, non_target_in_module = c, odds_ratio = unname(ft$estimate), p_value = ft$p.value)
}))
enrichment$BH_FDR_across_WGCNA_modules <- bh(enrichment$p_value)
write_csv(enrichment, file.path(OUTPUT, "tcga_coad_read_wgcna_target_module_enrichment.csv"))

module_summary <- aggregate(gene_symbol ~ module, data = gene_membership, FUN = length)
colnames(module_summary)[2] <- "module_size"
module_summary$target_overlap <- vapply(module_summary$module, function(m) sum(network_targets %in% gene_membership$gene_symbol[gene_membership$module == m]), integer(1))
module_summary$target_genes <- vapply(module_summary$module, function(m) paste(network_targets[network_targets %in% gene_membership$gene_symbol[gene_membership$module == m]], collapse = ";"), character(1))
module_summary <- merge(module_summary, module_trait, by = "module", all.x = TRUE, sort = FALSE)
write_csv(module_summary, file.path(OUTPUT, "tcga_coad_read_wgcna_module_summary.csv"))

pdf(file.path(OUTPUT, "tcga_coad_read_wgcna_module_trait_heatmap.pdf"), width = 12, height = 6)
cor_mat <- rbind(module_trait$target_score_cor, module_trait$cohort_COAD_cor)
rownames(cor_mat) <- c("41-gene target score", "COAD cohort")
colnames(cor_mat) <- module_trait$module
image(seq_len(ncol(cor_mat)), seq_len(nrow(cor_mat)), t(cor_mat), axes = FALSE, col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101), zlim = c(-1, 1), xlab = "WGCNA module", ylab = "Trait")
axis(1, at = seq_len(ncol(cor_mat)), labels = colnames(cor_mat), las = 2, cex.axis = 0.65)
axis(2, at = seq_len(nrow(cor_mat)), labels = rownames(cor_mat), las = 2)
box()
for (i in seq_len(nrow(cor_mat))) for (j in seq_len(ncol(cor_mat))) text(j, i, sprintf("%.2f", cor_mat[i, j]), cex = 0.7)
dev.off()

non_grey <- enrichment$module != "grey"
sig_modules <- enrichment[non_grey & enrichment$BH_FDR_across_WGCNA_modules < 0.05, , drop = FALSE]
top_trait <- module_trait[which.max(abs(module_trait$target_score_cor)), , drop = FALSE]
report <- c(
  "# TCGA-COAD + TCGA-READ full-transcriptome WGCNA", "",
  "## Scope", "",
  sprintf("A %s, %s WGCNA was built from the top %d genes selected by full-transcriptome MAD after gene-symbol deduplication. The fresh 41-gene DINP–CRC set was overlaid after network construction; target genes were not used to select or force network genes.", NETWORK_TYPE, COR_TYPE, ncol(datExpr)), "",
  sprintf("- Samples: **%d** primary tumors (%d COAD, %d READ).", nrow(datExpr), sum(metadata$cohort == "COAD"), sum(metadata$cohort == "READ")),
  sprintf("- Network genes: **%d**; fresh 41-gene targets present as network genes: **%d/41**.", ncol(datExpr), length(network_targets)),
  sprintf("- Targets projected to module eigengenes only: **%d/41**.", sum(!target_overlay$in_network_input)),
  sprintf("- Non-grey modules: **%d**; target-enriched modules by module-level BH-FDR<0.05: **%d**.", sum(non_grey), nrow(sig_modules)),
  sprintf("- Soft-threshold power: **%s** (%s).", soft_power, power_rule), "",
  "## Interpretation boundary", "",
  "WGCNA describes co-expression structure in CRC bulk tissue. It does not establish DINP exposure, direct chemical binding, causality, or direction of regulation. Module projections for targets absent from the variance-filtered network are descriptive and are not counted in Fisher enrichment. Bulk modules may also reflect tissue composition; macrophage claims require independent cell-type support.", "",
  "## Main outputs", "",
  "- `tcga_coad_read_wgcna_gene_module_membership.csv`: module labels and eigengene correlations for all network genes.",
  "- `tcga_coad_read_wgcna_41_gene_overlay.csv`: exact target overlay and projection-only results.",
  "- `tcga_coad_read_wgcna_target_module_enrichment.csv`: correct Fisher ORA using the observed network universe.",
  "- `tcga_coad_read_wgcna_module_trait_correlations.csv`: module eigengene associations with the 41-gene score and cohort.",
  "- `tcga_coad_read_wgcna_module_summary.csv`: module sizes, target overlap, and trait associations.",
  "- `tcga_coad_read_wgcna_sample_dendrogram.pdf`, `tcga_coad_read_wgcna_module_dendrogram.pdf`, and `tcga_coad_read_wgcna_module_trait_heatmap.pdf`: QC/overview figures.", "",
  "## Reproducibility", "",
  "- Expression scale: Xena-delivered log2(TPM+0.001).",
  sprintf("- Network: %s %s; power=%s; minModuleSize=30; mergeCutHeight=0.25; deepSplit=2.", NETWORK_TYPE, COR_TYPE, soft_power),
  sprintf("- R: %s; WGCNA: %s; dynamicTreeCut: %s; fastcluster: %s.", R.version.string, as.character(packageVersion("WGCNA")), as.character(packageVersion("dynamicTreeCut")), as.character(packageVersion("fastcluster"))),
  "- No automatic sample removal was performed; sample clustering is saved for inspection.",
  "- The 41-gene set was not used for network-gene selection."
)
writeLines(report, file.path(OUTPUT, "TCGA_COAD_READ_WGCNA_REPORT.md"), useBytes = TRUE)

manifest <- list(
  analysis = "TCGA COAD + READ full-transcriptome WGCNA with fresh 41-gene overlay",
  run_timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  input_matrix = normalizePath(INPUT_MATRIX, winslash = "/"), input_selection = normalizePath(INPUT_SELECTION, winslash = "/"), target_file = normalizePath(TARGET_FILE, winslash = "/"),
  samples = nrow(datExpr), coad_samples = sum(metadata$cohort == "COAD"), read_samples = sum(metadata$cohort == "READ"),
  network_genes = ncol(datExpr), target_genes = length(target_genes), target_genes_present_in_network = length(network_targets), target_genes_projection_only = sum(!target_overlay$in_network_input),
  module_count_non_grey = sum(non_grey), soft_threshold_power = soft_power, soft_threshold_rule = power_rule, target_enrichment_significant_modules = nrow(sig_modules),
  parameters = list(networkType = NETWORK_TYPE, corType = COR_TYPE, maxPOutliers = 0.05, minModuleSize = 30, mergeCutHeight = 0.25, deepSplit = 2),
  interpretation_boundary = "coexpression structure and post-network target overlay; not exposure causality or direct target confirmation",
  package_versions = list(R = R.version.string, WGCNA = as.character(packageVersion("WGCNA")), dynamicTreeCut = as.character(packageVersion("dynamicTreeCut")), fastcluster = as.character(packageVersion("fastcluster")))
)
jsonlite::write_json(manifest, file.path(OUTPUT, "tcga_coad_read_wgcna_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
message("WGCNA complete: ", nrow(datExpr), " samples; ", ncol(datExpr), " network genes; ", sum(non_grey), " non-grey modules")
