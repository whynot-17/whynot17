#!/usr/bin/env Rscript

# DINP-CRC project: reproduce transcriptomic preprocessing used in
# Zhu et al., Discover Oncology 2025 (DEHP-CRC), with DINP as the exposure.
#
# Outputs:
#   outputs/dinp_crc_dehp_replication/tcga_coread_expression_gene.tsv.gz
#   outputs/dinp_crc_dehp_replication/tcga_coread_pheno.tsv
#   outputs/dinp_crc_dehp_replication/GSE32323_gene_expression.tsv.gz
#   outputs/dinp_crc_dehp_replication/GSE32323_pheno.tsv
#   outputs/dinp_crc_dehp_replication/GSE21510_gene_expression.tsv.gz
#   outputs/dinp_crc_dehp_replication/GSE21510_pheno.tsv
#   outputs/dinp_crc_dehp_replication/GEO_merged_uncorrected.tsv.gz
#   outputs/dinp_crc_dehp_replication/GEO_merged_combat.tsv.gz
#   outputs/dinp_crc_dehp_replication/GEO_merged_pheno.tsv
#   outputs/dinp_crc_dehp_replication/preprocess_qc.tsv

options(stringsAsFactors = FALSE)
set.seed(20260908)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
outdir <- file.path(root, "outputs", "dinp_crc_dehp_replication")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

auto_bioc <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install(miss, ask = FALSE, update = FALSE)
  }
}
auto_cran <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
}

auto_bioc(c("GEOquery", "Biobase", "hgu133plus2.db", "AnnotationDbi", "sva"))
auto_cran(c("data.table"))

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(hgu133plus2.db)
  library(AnnotationDbi)
  library(sva)
  library(data.table)
})

write_tsv_gz <- function(x, path, rowname_col = NULL) {
  y <- as.data.frame(x, check.names = FALSE)
  if (!is.null(rowname_col)) y <- cbind(setNames(data.frame(rownames(y), check.names = FALSE), rowname_col), y)
  data.table::fwrite(y, file = path, sep = "\t", quote = FALSE)
}

collapse_probes_to_gene <- function(expr) {
  probe <- rownames(expr)
  map <- AnnotationDbi::select(hgu133plus2.db,
                               keys = probe,
                               keytype = "PROBEID",
                               columns = "SYMBOL")
  map <- map[!is.na(map$SYMBOL) & nzchar(map$SYMBOL), , drop = FALSE]
  map <- map[!duplicated(map$PROBEID), , drop = FALSE]
  idx <- match(map$PROBEID, rownames(expr))
  m <- expr[idx, , drop = FALSE]
  rownames(m) <- map$SYMBOL
  # reproducible gene-level collapse: median across probes mapping to same symbol
  genes <- unique(rownames(m))
  out <- vapply(genes, function(g) {
    z <- m[rownames(m) == g, , drop = FALSE]
    if (nrow(z) == 1L) as.numeric(z[1, ]) else apply(z, 2, median, na.rm = TRUE)
  }, numeric(ncol(m)))
  out <- t(out)
  colnames(out) <- colnames(expr)
  out
}

all_meta_text <- function(pd) {
  cols <- grep("title|source_name|characteristics", names(pd), ignore.case = TRUE, value = TRUE)
  apply(pd[, cols, drop = FALSE], 1, function(z) paste(z, collapse = " | "))
}

classify_crc <- function(pd, accession) {
  txt <- tolower(all_meta_text(pd))
  grp <- rep(NA_character_, length(txt))

  # normal first so 'non-cancerous' is never accidentally labeled tumor
  grp[grepl("non[- ]?cancer|non[- ]?tumou?r|adjacent normal|normal (colon|colorectal|mucosa|tissue)|normal tissue", txt)] <- "Normal"
  grp[is.na(grp) & grepl("colorectal cancer|colon cancer|rectal cancer|carcinoma|cancer tissue|tumou?r tissue|primary tumou?r|crc tissue", txt)] <- "Tumor"

  # GSE32323 contains additional 5-aza-treated/control cell-line samples: exclude explicitly.
  if (accession == "GSE32323") {
    is_cell <- grepl("cell line|5[- ]?aza|colo320|hct116|ht29|rko|sw480", txt)
    grp[is_cell] <- "Exclude_cell_line"
  }

  grp
}

process_geo <- function(acc) {
  message("Downloading ", acc, " ...")
  g <- GEOquery::getGEO(acc, GSEMatrix = TRUE, getGPL = FALSE)
  if (length(g) != 1L) stop(acc, ": expected one platform, got ", length(g))
  eset <- g[[1]]
  expr <- Biobase::exprs(eset)
  pd <- Biobase::pData(eset)
  pd$sample_id <- rownames(pd)
  pd$group <- classify_crc(pd, acc)

  keep <- pd$group %in% c("Normal", "Tumor")
  if (!any(keep)) stop(acc, ": group classifier found zero tissue samples. Inspect phenotype metadata.")
  expr <- expr[, pd$sample_id[keep], drop = FALSE]
  pd <- pd[keep, , drop = FALSE]

  # GEO records already report RMA/log2 processing for both GPL570 studies.
  q <- quantile(as.numeric(expr), probs = c(0, .25, .5, .75, 1), na.rm = TRUE)
  if (q[5] > 100) expr <- log2(expr + 1)

  gene <- collapse_probes_to_gene(expr)
  gene <- gene[rowSums(is.finite(gene)) == ncol(gene), , drop = FALSE]

  write_tsv_gz(gene, file.path(outdir, paste0(acc, "_gene_expression.tsv.gz")), "gene")
  data.table::fwrite(pd, file.path(outdir, paste0(acc, "_pheno.tsv")), sep = "\t", quote = FALSE)

  list(expr = gene, pheno = pd)
}

g32323 <- process_geo("GSE32323")
g21510 <- process_geo("GSE21510")

# Merge GEO cohorts at shared gene symbols.
common <- intersect(rownames(g32323$expr), rownames(g21510$expr))
x1 <- g32323$expr[common, , drop = FALSE]
x2 <- g21510$expr[common, , drop = FALSE]
geo_unc <- cbind(x1, x2)
geo_pheno <- rbind(
  data.frame(sample_id = colnames(x1), group = g32323$pheno[colnames(x1), "group"], batch = "GSE32323"),
  data.frame(sample_id = colnames(x2), group = g21510$pheno[colnames(x2), "group"], batch = "GSE21510")
)
rownames(geo_pheno) <- geo_pheno$sample_id

# Preserve an uncorrected merge (closest to the paper) and a ComBat version for robust downstream use.
write_tsv_gz(geo_unc, file.path(outdir, "GEO_merged_uncorrected.tsv.gz"), "gene")
data.table::fwrite(geo_pheno, file.path(outdir, "GEO_merged_pheno.tsv"), sep = "\t", quote = FALSE)

mod <- model.matrix(~ group, data = geo_pheno[colnames(geo_unc), , drop = FALSE])
geo_combat <- sva::ComBat(dat = geo_unc,
                          batch = geo_pheno[colnames(geo_unc), "batch"],
                          mod = mod,
                          par.prior = TRUE,
                          prior.plots = FALSE)
write_tsv_gz(geo_combat, file.path(outdir, "GEO_merged_combat.tsv.gz"), "gene")

# ---- TCGA COADREAD ----
# Reproduce the paper's TCGA-COADREAD concept through UCSC Xena legacy HiSeqV2.
# Dataset: TCGA.COADREAD.sampleMap/HiSeqV2
# We use direct hub downloads so sample barcodes determine tumor/normal status.
message("Downloading TCGA COADREAD HiSeqV2 from UCSC Xena ...")
tcga_url <- "https://tcga.xenahubs.net/download/TCGA.COADREAD.sampleMap/HiSeqV2.gz"
tcga_file <- tempfile(fileext = ".gz")
tryCatch(download.file(tcga_url, tcga_file, mode = "wb", quiet = TRUE),
         error = function(e) stop("TCGA Xena download failed: ", conditionMessage(e)))

tc <- data.table::fread(tcga_file, data.table = FALSE, check.names = FALSE)
if (!"sample" %in% names(tc)) names(tc)[1] <- "gene"
if (names(tc)[1] != "gene") names(tc)[1] <- "gene"
rownames(tc) <- tc$gene
tc$gene <- NULL
tc_mat <- as.matrix(tc)
storage.mode(tc_mat) <- "numeric"

# TCGA barcode sample-type code positions 14-15: 01 primary tumor, 11 solid tissue normal.
barcodes <- colnames(tc_mat)
type_code <- substr(barcodes, 14, 15)
keep <- type_code %in% c("01", "11")
tc_mat <- tc_mat[, keep, drop = FALSE]
type_code <- type_code[keep]
tc_pheno <- data.frame(
  sample_id = colnames(tc_mat),
  group = ifelse(type_code == "11", "Normal", "Tumor"),
  sample_type_code = type_code,
  stringsAsFactors = FALSE
)
rownames(tc_pheno) <- tc_pheno$sample_id

# HiSeqV2 is already normalized/log transformed; average duplicate gene symbols if any.
valid <- !is.na(rownames(tc_mat)) & nzchar(rownames(tc_mat))
tc_mat <- tc_mat[valid, , drop = FALSE]
if (anyDuplicated(rownames(tc_mat))) {
  genes <- unique(rownames(tc_mat))
  tmp <- vapply(genes, function(g) {
    z <- tc_mat[rownames(tc_mat) == g, , drop = FALSE]
    if (nrow(z) == 1L) as.numeric(z[1, ]) else colMeans(z, na.rm = TRUE)
  }, numeric(ncol(tc_mat)))
  tc_mat <- t(tmp)
  colnames(tc_mat) <- tc_pheno$sample_id
}

write_tsv_gz(tc_mat, file.path(outdir, "tcga_coread_expression_gene.tsv.gz"), "gene")
data.table::fwrite(tc_pheno, file.path(outdir, "tcga_coread_pheno.tsv"), sep = "\t", quote = FALSE)

qc <- rbind(
  data.frame(dataset = "GSE32323", Normal = sum(g32323$pheno$group == "Normal"), Tumor = sum(g32323$pheno$group == "Tumor"), Genes = nrow(g32323$expr)),
  data.frame(dataset = "GSE21510", Normal = sum(g21510$pheno$group == "Normal"), Tumor = sum(g21510$pheno$group == "Tumor"), Genes = nrow(g21510$expr)),
  data.frame(dataset = "GEO_merged", Normal = sum(geo_pheno$group == "Normal"), Tumor = sum(geo_pheno$group == "Tumor"), Genes = nrow(geo_unc)),
  data.frame(dataset = "TCGA_COADREAD", Normal = sum(tc_pheno$group == "Normal"), Tumor = sum(tc_pheno$group == "Tumor"), Genes = nrow(tc_mat))
)
data.table::fwrite(qc, file.path(outdir, "preprocess_qc.tsv"), sep = "\t", quote = FALSE)
print(qc)

expected <- data.frame(dataset = c("GSE32323", "GSE21510", "GEO_merged", "TCGA_COADREAD"),
                       Normal = c(17, 25, 42, 41),
                       Tumor = c(17, 123, 140, 476))
message("Paper-reported target counts:")
print(expected)
message("NOTE: if current TCGA Xena counts differ from 41/476, do not silently subset to force agreement; audit the paper's exact TCGA release/sample filtering first.")
