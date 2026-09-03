# DINP--CRC bulk co-expression analysis (standardized v2)
#
# Primary purpose: build an unbiased CRC bulk co-expression network and then
# overlay the frozen 81-gene DINP--CRC intersection and the 7-gene macrophage
# driver set.  The script uses
# the public, processed GSE39582 GPL570 series matrix and keeps two analyses:
# (1) all samples, including the 19 non-tumoral samples, for tumor-status
# module association; and (2) tumor-only samples, to avoid making the network
# entirely a tumor-versus-normal contrast.
#
# The script deliberately does not infer causality from co-expression.  It
# reports module membership, post hoc target-set overlay/enrichment, a fixed
# macrophage marker-score association, and module-trait associations as
# exploratory convergence evidence. Frozen target/driver genes are never
# forced into the WGCNA input.

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  hit <- which(args == flag)
  if (length(hit) == 0L || hit[1] == length(args)) return(default)
  args[hit[1] + 1L]
}

matrix_path <- get_arg("--matrix", "D:/whynot17/work/wgcna_crc/raw/GSE39582_series_matrix.txt.gz")
annotation_path <- get_arg("--annotation", "D:/whynot17/work/wgcna_crc/raw/GPL570.annot.gz")
target_path <- get_arg(
  "--targets",
  "analysis/dinp_crc_81gene_singlecell_localization/outputs/input_81_genes.csv"
)
driver_path <- get_arg(
  "--drivers",
  "analysis/dinp_crc_81gene_macrophage_driver_decomposition/outputs/macrophage_driver_candidates.csv"
)
out_dir <- get_arg("--outdir", "analysis/dinp_crc_wgcna/outputs")
max_genes <- as.integer(get_arg("--max-genes", "5000"))
seed <- as.integer(get_arg("--seed", "39582"))

# A fixed, canonical macrophage/myeloid marker panel is used only to create a
# bulk-expression abundance/state proxy. It is not a deconvolution result.
macrophage_markers <- c(
  "CD68", "LST1", "TYROBP", "AIF1", "FCER1G", "CTSS", "C1QA", "C1QB",
  "C1QC", "MS4A7", "LILRB1", "LGALS3", "CD14", "CTSB", "SPI1"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# The WGCNA installation is intentionally kept off C: with the project data.
project_rlib <- "D:/whynot17/Rlib"
if (dir.exists(project_rlib)) .libPaths(c(project_rlib, .libPaths()))
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  stop("WGCNA is not available. Expected it in D:/whynot17/Rlib.")
}
suppressPackageStartupMessages(library(WGCNA))

if (!file.exists(matrix_path)) stop("Missing expression matrix: ", matrix_path)
if (!file.exists(annotation_path)) stop("Missing GPL570 annotation: ", annotation_path)
if (!file.exists(target_path)) stop("Missing 81-gene list: ", target_path)
if (!file.exists(driver_path)) stop("Missing macrophage driver list: ", driver_path)

strip_outer_quotes <- function(x) {
  x <- as.character(x)
  x <- sub('^"', "", x)
  x <- sub('"$', "", x)
  x
}

geo_tab_values <- function(line) {
  parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
  if (length(parts) <= 1L) return(character())
  strip_outer_quotes(parts[-1L])
}

read_geo_metadata <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  lines <- character()
  repeat {
    ln <- readLines(con, n = 1L, warn = FALSE)
    if (length(ln) == 0L) break
    if (grepl("!series_matrix_table_begin", ln, fixed = TRUE)) break
    lines <- c(lines, ln)
  }

  accession_line <- lines[startsWith(lines, "!Sample_geo_accession")][1]
  source_line <- lines[startsWith(lines, "!Sample_source_name_ch1")][1]
  if (is.na(accession_line) || is.na(source_line)) {
    stop("GEO matrix metadata is missing sample accession/source lines.")
  }
  sample_id <- geo_tab_values(accession_line)
  source_name <- geo_tab_values(source_line)
  meta <- data.frame(
    sample_id = sample_id,
    source_name = source_name,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  characteristic_lines <- lines[startsWith(lines, "!Sample_characteristics_ch1")]
  if (length(characteristic_lines) > 0L) {
    for (ln in characteristic_lines) {
      vals <- geo_tab_values(ln)
      if (length(vals) != nrow(meta)) next
      first_nonempty <- vals[nzchar(vals)][1]
      if (is.na(first_nonempty) || !grepl(":", first_nonempty, fixed = TRUE)) next
      key <- trimws(sub(":.*$", "", first_nonempty))
      if (!nzchar(key)) next
      value <- trimws(sub("^[^:]*:", "", vals))
      # Duplicate keys are uncommon in this series.  If present, preserve the
      # first and make later columns explicit rather than silently overwriting.
      col_name <- key
      if (col_name %in% names(meta)) {
        suffix <- 2L
        while (paste0(col_name, "_", suffix) %in% names(meta)) suffix <- suffix + 1L
        col_name <- paste0(col_name, "_", suffix)
      }
      meta[[col_name]] <- value
    }
  }

  meta$is_non_tumor <- grepl("non[- ]?tumor|non[- ]?tumoral|normal", tolower(meta$source_name))
  meta$tumor_status <- ifelse(meta$is_non_tumor, "non_tumoral", "tumor")
  meta
}

read_geo_expression <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  repeat {
    ln <- readLines(con, n = 1L, warn = FALSE)
    if (length(ln) == 0L) stop("GEO matrix table begin marker not found.")
    if (grepl("!series_matrix_table_begin", ln, fixed = TRUE)) break
  }
  tbl <- read.delim(
    con,
    header = TRUE,
    sep = "\t",
    quote = "\"",
    comment.char = "",
    fill = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (ncol(tbl) < 3L) stop("Expression table has unexpectedly few columns.")
  first_col <- names(tbl)[1]
  keep <- !grepl("^!series_matrix_table_end", as.character(tbl[[1]]))
  tbl <- tbl[keep, , drop = FALSE]
  probe <- as.character(tbl[[1]])
  values <- as.matrix(tbl[, -1L, drop = FALSE])
  suppressWarnings(storage.mode(values) <- "numeric")
  rownames(values) <- probe
  colnames(values) <- names(tbl)[-1L]
  list(matrix = values, sample_ids = colnames(values), first_col = first_col)
}

read_platform_annotation <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  repeat {
    ln <- readLines(con, n = 1L, warn = FALSE)
    if (length(ln) == 0L) stop("GPL annotation table begin marker not found.")
    if (grepl("!platform_table_begin", ln, fixed = TRUE)) break
  }
  ann <- read.delim(
    con,
    header = TRUE,
    sep = "\t",
    quote = "\"",
    comment.char = "",
    fill = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (!all(c("ID", "Gene symbol") %in% names(ann))) {
    stop("GPL570 annotation lacks expected ID/Gene symbol columns.")
  }
  ann <- ann[!grepl("^!platform_table_end", as.character(ann$ID)), , drop = FALSE]
  ann
}

message("Reading GEO metadata ...")
meta <- read_geo_metadata(matrix_path)
message("Reading GEO expression matrix ...")
geo <- read_geo_expression(matrix_path)
expr_probe <- geo$matrix

if (!all(meta$sample_id %in% colnames(expr_probe))) {
  missing_meta <- setdiff(meta$sample_id, colnames(expr_probe))
  stop("Metadata samples missing from expression matrix: ", paste(head(missing_meta, 10), collapse = ", "))
}
meta <- meta[match(colnames(expr_probe), meta$sample_id), , drop = FALSE]
rownames(meta) <- meta$sample_id

message("Reading GPL570 probe annotation ...")
ann <- read_platform_annotation(annotation_path)
ann$probe_id <- as.character(ann$ID)
ann$gene_symbol_raw <- trimws(as.character(ann[["Gene symbol"]]))
ann$gene_symbol <- vapply(strsplit(ann$gene_symbol_raw, "///", fixed = TRUE), function(x) {
  x <- trimws(x)
  x <- x[nzchar(x)]
  if (length(x) == 0L) return(NA_character_)
  x[1L]
}, character(1))
ann$gene_symbol[ann$gene_symbol %in% c("", "---", "NA", "N/A")] <- NA_character_
ann <- ann[match(rownames(expr_probe), ann$probe_id), , drop = FALSE]
ann$probe_id <- rownames(expr_probe)

target_tbl <- read.csv(target_path, check.names = FALSE, stringsAsFactors = FALSE)
targets <- unique(trimws(as.character(target_tbl$gene_symbol)))
targets <- targets[nzchar(targets)]
driver_tbl <- read.csv(driver_path, check.names = FALSE, stringsAsFactors = FALSE)
drivers <- unique(trimws(as.character(driver_tbl$gene_symbol)))
drivers <- drivers[nzchar(drivers)]

valid_mapping <- !is.na(ann$gene_symbol) & nzchar(ann$gene_symbol)
mapped_expr <- expr_probe[valid_mapping, , drop = FALSE]
mapped_gene <- ann$gene_symbol[valid_mapping]
mapped_probe <- rownames(mapped_expr)

gene_groups <- split(seq_len(nrow(mapped_expr)), mapped_gene)
gene_ids <- names(gene_groups)
gene_expr <- matrix(
  NA_real_,
  nrow = length(gene_groups),
  ncol = ncol(mapped_expr),
  dimnames = list(gene_ids, colnames(mapped_expr))
)
selected_probe <- character(length(gene_groups))
names(selected_probe) <- gene_ids
for (i in seq_along(gene_groups)) {
  idx <- gene_groups[[i]]
  if (length(idx) == 1L) {
    pick <- idx
  } else {
    probe_mad <- apply(mapped_expr[idx, , drop = FALSE], 1L, mad, na.rm = TRUE)
    probe_mad[is.na(probe_mad)] <- -Inf
    pick <- idx[which.max(probe_mad)]
  }
  gene_expr[i, ] <- as.numeric(mapped_expr[pick, ])
  selected_probe[i] <- mapped_probe[pick]
}

# GEO processed matrices should be complete, but median-impute any residual
# missing values before WGCNA and record the event in the audit.
missing_before_impute <- sum(is.na(gene_expr))
if (missing_before_impute > 0L) {
  for (i in seq_len(nrow(gene_expr))) {
    miss <- is.na(gene_expr[i, ])
    if (any(miss)) {
      med <- median(gene_expr[i, !miss], na.rm = TRUE)
      gene_expr[i, miss] <- med
    }
  }
}

gene_mad <- apply(gene_expr, 1L, mad, na.rm = TRUE)
gene_mad[is.na(gene_mad)] <- 0
ranked_genes <- names(sort(gene_mad, decreasing = TRUE))
target_present <- intersect(targets, rownames(gene_expr))
driver_present <- intersect(drivers, rownames(gene_expr))
selected_genes <- head(ranked_genes, max_genes)
selected_genes <- selected_genes[selected_genes %in% rownames(gene_expr)]
target_selected <- intersect(targets, selected_genes)
driver_selected <- intersect(drivers, selected_genes)
macrophage_present <- intersect(macrophage_markers, rownames(gene_expr))

write.csv(
  data.frame(
    gene_symbol = targets,
    present_in_gpl570 = targets %in% rownames(gene_expr),
    selected_by_mad_filter = targets %in% selected_genes,
    selected_probe = unname(selected_probe[targets]),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "target_gene_mapping_audit.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  data.frame(
    gene_symbol = drivers,
    present_in_gpl570 = drivers %in% rownames(gene_expr),
    selected_by_mad_filter = drivers %in% selected_genes,
    selected_probe = unname(selected_probe[drivers]),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "macrophage_driver_mapping_audit.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  data.frame(
    gene_symbol = macrophage_markers,
    present_in_gene_matrix = macrophage_markers %in% rownames(gene_expr),
    selected_by_mad_filter = macrophage_markers %in% selected_genes,
    selected_probe = unname(selected_probe[macrophage_markers]),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "macrophage_marker_mapping_audit.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  data.frame(
    marker_set = "fixed_macrophage_myeloid_proxy",
    gene_symbol = macrophage_markers,
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "macrophage_core_marker_set.csv"),
  row.names = FALSE
)
write.csv(data.frame(gene_symbol = selected_genes), file.path(out_dir, "wgcna_selected_genes.csv"), row.names = FALSE)

parse_numeric <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "null", "NULL")] <- NA_character_
  suppressWarnings(as.numeric(gsub("[^0-9.-]", "", x)))
}

parse_stage <- function(x) {
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA_real_, length(y))
  out[grepl("stage\\s*(iv|4)|\\biv\\b", y)] <- 4
  out[grepl("stage\\s*(iii|3)|\\biii\\b", y)] <- 3
  out[grepl("stage\\s*(ii|2)|\\bii\\b", y)] <- 2
  out[grepl("stage\\s*(i|1)|\\bi\\b", y)] <- 1
  out
}

add_factor_dummies <- function(traits, meta_sub, key, prefix, min_n = 5L) {
  if (!key %in% names(meta_sub)) return(traits)
  x <- trimws(as.character(meta_sub[[key]]))
  x[x %in% c("", "NA", "N/A", "null", "NULL")] <- NA_character_
  tab <- sort(table(x), decreasing = TRUE)
  keep <- names(tab)[tab >= min_n]
  if (length(keep) < 2L) return(traits)
  # A reference category is omitted to prevent a redundant all-category set.
  keep <- keep[-1L]
  for (lev in keep) traits[[paste0(prefix, "_", make.names(lev))]] <- as.numeric(x == lev)
  traits
}

build_traits <- function(meta_sub, include_status = TRUE, macrophage_score = NULL) {
  traits <- data.frame(row.names = meta_sub$sample_id, check.names = FALSE)
  if (include_status) traits$tumor_status <- as.numeric(meta_sub$tumor_status == "tumor")
  if ("age.at.diagnosis (year)" %in% names(meta_sub)) {
    traits$age_at_diagnosis <- parse_numeric(meta_sub[["age.at.diagnosis (year)"]])
  }
  if ("tnm.stage" %in% names(meta_sub)) traits$stage_ordinal <- parse_stage(meta_sub$tnm.stage)
  if ("os.event" %in% names(meta_sub)) traits$os_event <- parse_numeric(meta_sub$os.event)
  if ("rfs.event" %in% names(meta_sub)) traits$rfs_event <- parse_numeric(meta_sub$rfs.event)
  for (spec in list(
    c("dataset", "dataset"),
    c("Sex", "sex"),
    c("tumor.location", "tumor_location"),
    c("cit.molecularsubtype", "molecular_subtype"),
    c("mmr.status", "mmr_status"),
    c("cimp.status", "cimp_status"),
    c("cin.status", "cin_status")
  )) traits <- add_factor_dummies(traits, meta_sub, spec[1], spec[2])
  if (!is.null(macrophage_score)) {
    traits$macrophage_core_marker_score <- unname(macrophage_score[meta_sub$sample_id])
  }
  traits <- traits[, vapply(traits, function(x) sum(!is.na(x)) >= 10L && length(unique(x[!is.na(x)])) >= 2L, logical(1)), drop = FALSE]
  traits
}

cor_pvalue <- function(r, n) {
  if (is.na(r) || n < 4L || abs(r) >= 1) return(ifelse(is.na(r), NA_real_, 0))
  WGCNA::corPvalueStudent(r, n)
}

run_one_network <- function(label, sample_idx, include_status) {
  message("Preparing ", label, " network ...")
  gene_mat <- gene_expr[selected_genes, sample_idx, drop = FALSE]
  datExpr <- t(gene_mat)
  colnames(datExpr) <- selected_genes
  rownames(datExpr) <- colnames(gene_mat)

  gsg <- goodSamplesGenes(datExpr, verbose = 0)
  if (!gsg$allOK) {
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  }
  if (nrow(datExpr) < 50L || ncol(datExpr) < 100L) {
    stop(label, " network has too few samples or genes after QC: ", nrow(datExpr), " x ", ncol(datExpr))
  }

  powers <- c(1:10, seq(12, 30, by = 2))
  message("Selecting soft threshold for ", label, " ...")
  sft <- suppressWarnings(pickSoftThreshold(
    datExpr,
    powerVector = powers,
    networkType = "signed",
    corFnc = "bicor",
    corOptions = list(maxPOutliers = 0.05),
    verbose = 0
  ))
  fit <- sft$fitIndices
  r2 <- fit$SFT.R.sq
  candidate <- which(!is.na(r2) & r2 >= 0.80)
  if (length(candidate) > 0L) {
    soft_power <- fit$Power[min(candidate)]
    power_rule <- "smallest power with scale-free fit R2 >= 0.80"
  } else {
    # Transparent fallback: maximize the observed fit, with lower power on ties.
    best <- which.max(ifelse(is.na(r2), -Inf, r2))
    soft_power <- fit$Power[best]
    power_rule <- "maximum observed scale-free fit; no power reached R2 >= 0.80"
  }
  fit$selected <- fit$Power == soft_power
  fit$analysis <- label
  write.csv(fit, file.path(out_dir, paste0("soft_threshold_", label, ".csv")), row.names = FALSE)

  message("Building ", label, " modules at power ", soft_power, " ...")
  set.seed(seed)
  net <- blockwiseModules(
    datExpr,
    power = soft_power,
    networkType = "signed",
    TOMType = "signed",
    minModuleSize = 30,
    reassignThreshold = 0,
    mergeCutHeight = 0.25,
    deepSplit = 2,
    pamRespectsDendro = FALSE,
    numericLabels = TRUE,
    saveTOMs = FALSE,
    verbose = 2,
    corType = "bicor",
    maxPOutliers = 0.05,
    maxBlockSize = max(6000L, ncol(datExpr) + 10L)
  )
  colors <- labels2colors(net$colors)
  names(colors) <- colnames(datExpr)
  MEs <- orderMEs(net$MEs)
  # blockwiseModules(numericLabels=TRUE) names eigengenes ME0, ME1, ...,
  # whereas the gene assignments below are represented by WGCNA colors.
  # Translate the numeric eigengene labels before computing own-module kME.
  me_numeric <- suppressWarnings(as.integer(sub("^ME", "", colnames(MEs))))
  me_colors <- labels2colors(me_numeric)
  colnames(MEs) <- paste0("ME", me_colors)
  if (ncol(MEs) == 0L) stop("No non-grey modules were returned for ", label)
  MEs_non_grey <- MEs[, colnames(MEs) != "MEgrey", drop = FALSE]
  if (ncol(MEs_non_grey) == 0L) stop("No non-grey modules were returned for ", label)

  meta_sub <- meta[sample_idx, , drop = FALSE]
  rownames(meta_sub) <- meta_sub$sample_id
  macrophage_score <- NULL
  if (length(macrophage_present) >= 3L) {
    marker_mat <- t(gene_expr[macrophage_present, sample_idx, drop = FALSE])
    marker_z <- scale(marker_mat)
    macrophage_score <- rowMeans(marker_z, na.rm = TRUE)
    macrophage_score[!is.finite(macrophage_score)] <- NA_real_
    names(macrophage_score) <- rownames(marker_mat)
  }
  traits <- build_traits(meta_sub, include_status = include_status, macrophage_score = macrophage_score)
  traits <- traits[rownames(datExpr), , drop = FALSE]

  module_trait_cor <- matrix(NA_real_, nrow = ncol(MEs_non_grey), ncol = ncol(traits), dimnames = list(colnames(MEs_non_grey), colnames(traits)))
  module_trait_p <- module_trait_cor
  module_trait_n <- module_trait_cor
  if (ncol(traits) > 0L) {
    for (i in seq_len(ncol(MEs_non_grey))) {
      for (j in seq_len(ncol(traits))) {
        keep <- complete.cases(MEs_non_grey[, i], traits[, j])
        n <- sum(keep)
        r <- if (n >= 4L) suppressWarnings(cor(MEs_non_grey[keep, i], traits[keep, j], method = "pearson")) else NA_real_
        module_trait_cor[i, j] <- r
        module_trait_p[i, j] <- cor_pvalue(r, n)
        module_trait_n[i, j] <- n
      }
    }
  }
  write.csv(as.data.frame(module_trait_cor), file.path(out_dir, paste0("module_trait_correlation_", label, ".csv")))
  write.csv(as.data.frame(module_trait_p), file.path(out_dir, paste0("module_trait_pvalue_", label, ".csv")))
  write.csv(as.data.frame(module_trait_n), file.path(out_dir, paste0("module_trait_complete_n_", label, ".csv")))
  module_trait_fdr <- module_trait_p
  if (ncol(module_trait_fdr) > 0L) {
    for (j in seq_len(ncol(module_trait_fdr))) {
      module_trait_fdr[, j] <- p.adjust(module_trait_fdr[, j], method = "BH")
    }
  }
  write.csv(as.data.frame(module_trait_fdr), file.path(out_dir, paste0("module_trait_fdr_", label, ".csv")))

  if ("macrophage_core_marker_score" %in% colnames(module_trait_cor)) {
    mac_idx <- which(colnames(module_trait_cor) == "macrophage_core_marker_score")[1L]
    macrophage_assoc <- data.frame(
      analysis = label,
      module = rownames(module_trait_cor),
      macrophage_score_cor = unname(module_trait_cor[, mac_idx]),
      macrophage_score_p = unname(module_trait_p[, mac_idx]),
      macrophage_score_fdr = unname(module_trait_fdr[, mac_idx]),
      stringsAsFactors = FALSE
    )
  } else {
    macrophage_assoc <- data.frame(
      analysis = label,
      module = rownames(module_trait_cor),
      macrophage_score_cor = NA_real_,
      macrophage_score_p = NA_real_,
      macrophage_score_fdr = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  write.csv(macrophage_assoc, file.path(out_dir, paste0("macrophage_module_association_", label, ".csv")), row.names = FALSE, na = "")

  if (ncol(traits) > 0L) {
    png(file.path(out_dir, paste0("module_trait_heatmap_", label, ".png")), width = 2200, height = max(1600, 80 * nrow(module_trait_cor)), res = 220)
    labeledHeatmap(
      Matrix = module_trait_cor,
      xLabels = colnames(module_trait_cor),
      yLabels = rownames(module_trait_cor),
      ySymbols = rownames(module_trait_cor),
      colorLabels = FALSE,
      colors = blueWhiteRed(50),
      textMatrix = ifelse(is.na(module_trait_p), "", ifelse(module_trait_p < 0.001, "***", ifelse(module_trait_p < 0.01, "**", ifelse(module_trait_p < 0.05, "*", "")))),
      setStdMargins = FALSE,
      cex.text = 0.55,
      zlim = c(-1, 1),
      main = paste("GSE39582", label, "module--trait correlations")
    )
    dev.off()
  }

  kME <- suppressWarnings(cor(datExpr, MEs, use = "pairwise.complete.obs", method = "pearson"))
  colnames(kME) <- sub("^ME", "kME_", colnames(kME))
  own_kme <- rep(NA_real_, nrow(kME))
  for (i in seq_len(nrow(kME))) {
    own_col <- paste0("kME_", colors[i])
    if (own_col %in% colnames(kME)) own_kme[i] <- kME[i, own_col]
  }
  membership <- data.frame(
    gene_symbol = rownames(kME),
    module = unname(colors),
    kME_own = own_kme,
    stringsAsFactors = FALSE
  )
  membership <- cbind(membership, as.data.frame(kME, check.names = FALSE))
  membership <- membership[order(membership$module, -abs(membership$kME_own)), , drop = FALSE]
  write.csv(membership, file.path(out_dir, paste0("gene_module_membership_", label, ".csv")), row.names = FALSE)

  target_map <- merge(
    data.frame(gene_symbol = targets, stringsAsFactors = FALSE),
    membership[, c("gene_symbol", "module", "kME_own"), drop = FALSE],
    by = "gene_symbol", all.x = TRUE, sort = FALSE
  )
  target_map <- target_map[match(targets, target_map$gene_symbol), , drop = FALSE]
  target_map$selected_in_wgcna_input <- target_map$gene_symbol %in% colnames(datExpr)
  target_map$target_group <- ifelse(target_map$gene_symbol %in% drivers, "81_gene_program_and_7_driver", "81_gene_program")
  target_map <- target_map[order(is.na(target_map$module), target_map$module, -abs(target_map$kME_own)), , drop = FALSE]
  write.csv(target_map, file.path(out_dir, paste0("target_module_mapping_", label, ".csv")), row.names = FALSE)

  # Enrichment is evaluated only for target/driver genes that were selected by
  # the prespecified variance filter. All frozen genes are retained in the
  # overlay table above, including genes absent from the natural WGCNA input.
  target_input <- intersect(targets, colnames(datExpr))
  driver_input <- intersect(drivers, colnames(datExpr))

  module_levels <- sort(unique(colors))
  enrich_rows <- lapply(module_levels, function(mod) {
    in_mod <- names(colors)[colors == mod]
    a <- sum(target_input %in% in_mod)
    b <- length(target_input) - a
    target_nonmod <- sum(!(names(colors) %in% target_input) & colors == mod)
    target_other <- sum(!(names(colors) %in% target_input) & colors != mod)
    mat <- matrix(c(a, b, target_nonmod, target_other), nrow = 2L, byrow = TRUE)
    p <- tryCatch(fisher.test(mat, alternative = "greater")$p.value, error = function(e) NA_real_)
    da <- sum(driver_input %in% in_mod)
    db <- length(driver_input) - da
    driver_nonmod <- sum(!(names(colors) %in% driver_input) & colors == mod)
    driver_other <- sum(!(names(colors) %in% driver_input) & colors != mod)
    dmat <- matrix(c(da, db, driver_nonmod, driver_other), nrow = 2L, byrow = TRUE)
    dp <- tryCatch(fisher.test(dmat, alternative = "greater")$p.value, error = function(e) NA_real_)
    data.frame(
      analysis = label,
      module = mod,
      module_gene_n = length(in_mod),
      target_gene_n_in_input = length(target_input),
      target_gene_n_in_module = a,
      target_enrichment_OR = tryCatch(unname(fisher.test(mat)$estimate), error = function(e) NA_real_),
      target_enrichment_p = p,
      target_enrichment_fdr = NA_real_,
      driver_gene_n_in_input = length(driver_input),
      driver_gene_n_in_module = da,
      driver_enrichment_fdr = NA_real_,
      driver_enrichment_p = dp,
      stringsAsFactors = FALSE
    )
  })
  enrichment <- do.call(rbind, enrich_rows)
  # Grey is WGCNA's unassigned bin, not a biological module.  Keep it in the
  # audit table but do not let it enter the module-level multiplicity family.
  enrichment$target_enrichment_fdr <- NA_real_
  enrichment$driver_enrichment_fdr <- NA_real_
  non_grey <- enrichment$module != "grey"
  enrichment$target_enrichment_fdr[non_grey] <- p.adjust(enrichment$target_enrichment_p[non_grey], method = "BH")
  enrichment$driver_enrichment_fdr[non_grey] <- p.adjust(enrichment$driver_enrichment_p[non_grey], method = "BH")
  enrichment <- enrichment[order(enrichment$target_enrichment_p, enrichment$module), , drop = FALSE]
  write.csv(enrichment, file.path(out_dir, paste0("module_target_enrichment_", label, ".csv")), row.names = FALSE)

  module_summary <- aggregate(gene_symbol ~ module, membership, length)
  names(module_summary)[2] <- "module_gene_n"
  target_map_selected <- target_map[!is.na(target_map$module), , drop = FALSE]
  target_counts <- aggregate(gene_symbol ~ module, target_map_selected, length)
  names(target_counts)[2] <- "target_gene_n"
  module_summary <- merge(module_summary, target_counts, by = "module", all.x = TRUE)
  module_summary$target_gene_n[is.na(module_summary$target_gene_n)] <- 0L
  module_summary <- merge(module_summary, enrichment[, c("module", "target_enrichment_p", "target_enrichment_fdr", "driver_gene_n_in_module", "driver_enrichment_p", "driver_enrichment_fdr")], by = "module", all.x = TRUE)
  module_summary <- module_summary[order(module_summary$target_enrichment_p, module_summary$module), , drop = FALSE]
  write.csv(module_summary, file.path(out_dir, paste0("module_summary_", label, ".csv")), row.names = FALSE)

  top_target <- target_map_selected[order(-abs(target_map_selected$kME_own)), , drop = FALSE]
  if (nrow(top_target) > 0L) write.csv(head(top_target, 30L), file.path(out_dir, paste0("top_target_membership_", label, ".csv")), row.names = FALSE)

  audit <- data.frame(
    analysis = label,
    input_sample_n = length(sample_idx),
    wgcna_sample_n = nrow(datExpr),
    wgcna_gene_n = ncol(datExpr),
    non_grey_module_n = ncol(MEs_non_grey),
    grey_gene_n = sum(colors == "grey"),
    soft_power = soft_power,
    soft_power_rule = power_rule,
    target_input_n = length(target_input),
    driver_input_n = length(driver_input),
    macrophage_marker_n = length(macrophage_present),
    macrophage_score_available = "macrophage_core_marker_score" %in% colnames(traits),
    stringsAsFactors = FALSE
  )
  write.csv(audit, file.path(out_dir, paste0("network_audit_", label, ".csv")), row.names = FALSE)

  list(
    label = label,
    n_samples = nrow(datExpr),
    n_genes = ncol(datExpr),
    n_modules = ncol(MEs_non_grey),
    soft_power = soft_power,
    power_rule = power_rule,
    module_trait_cor = module_trait_cor,
    module_trait_p = module_trait_p,
    enrichment = enrichment,
    target_map = target_map,
    macrophage_assoc = macrophage_assoc
  )
}

all_idx <- seq_len(ncol(gene_expr))
tumor_idx <- which(!meta$is_non_tumor)
if (length(tumor_idx) < 50L) stop("Tumor-only subset has too few samples: ", length(tumor_idx))

sample_audit <- data.frame(
  sample_n = nrow(meta),
  tumor_n = sum(!meta$is_non_tumor),
  non_tumoral_n = sum(meta$is_non_tumor),
  expression_probe_n = nrow(expr_probe),
  mapped_gene_n = nrow(gene_expr),
  selected_gene_n = length(selected_genes),
  target_n = length(targets),
  target_present_n = length(target_present),
  target_selected_n = length(target_selected),
  driver_n = length(drivers),
  driver_present_n = length(driver_present),
  driver_selected_n = length(driver_selected),
  macrophage_marker_n = length(macrophage_present),
  missing_expression_values_before_impute = missing_before_impute,
  source_matrix = "GSE39582_series_matrix.txt.gz",
  platform_annotation = "GPL570.annot.gz",
  stringsAsFactors = FALSE
)
write.csv(sample_audit, file.path(out_dir, "input_sample_gene_audit.csv"), row.names = FALSE)
write.csv(meta, file.path(out_dir, "sample_metadata_audit.csv"), row.names = FALSE)

all_result <- run_one_network("all_samples", all_idx, include_status = TRUE)
tumor_result <- run_one_network("tumor_only", tumor_idx, include_status = FALSE)

all_top <- head(all_result$enrichment[all_result$enrichment$module != "grey", , drop = FALSE], 5L)
tumor_top <- head(tumor_result$enrichment[tumor_result$enrichment$module != "grey", , drop = FALSE], 5L)
fmt <- function(x) ifelse(is.na(x), "NA", formatC(x, digits = 4, format = "g"))
report <- c(
  "# DINP--CRC bulk WGCNA report",
  "",
  "## Scope",
  "",
  "This analysis builds an unbiased CRC co-expression network and overlays the frozen 81-gene DINP--CRC intersection and seven macrophage driver candidates. WGCNA is used as exploratory co-expression convergence evidence; it does not establish that DINP causes any module or trait.",
  "",
  "## Input audit",
  "",
  paste0("- Samples: ", nrow(meta), " total; ", sum(!meta$is_non_tumor), " tumor and ", sum(meta$is_non_tumor), " non-tumoral."),
  paste0("- Probe sets: ", nrow(expr_probe), "; collapsed gene-level matrix: ", nrow(gene_expr), " genes."),
  paste0("- WGCNA input: ", length(selected_genes), " genes selected by MAD only; frozen target/driver genes were not forced into network construction."),
  paste0("- Frozen 81-gene list present on GPL570: ", length(target_present), "/", length(targets), "; selected by MAD: ", length(target_selected), "; seven-driver list present: ", length(driver_present), "/", length(drivers), "; selected by MAD: ", length(driver_selected), "."),
  paste0("- Macrophage/myeloid proxy: ", length(macrophage_present), "/", length(macrophage_markers), " fixed markers available; the resulting score is an abundance/state proxy, not cell deconvolution."),
  "- Duplicate probes were collapsed per gene by retaining the probe with the largest across-sample MAD.",
  "- The processed GEO matrix was used as supplied; no outcome-driven gene filtering was applied.",
  "",
  "## Network analyses",
  "",
  paste0("- All samples: n=", all_result$n_samples, "; ", all_result$n_modules, " non-grey modules; soft power=", all_result$soft_power, "."),
  paste0("- Tumor-only: n=", tumor_result$n_samples, "; ", tumor_result$n_modules, " non-grey modules; soft power=", tumor_result$soft_power, "."),
  "- The all-sample analysis enables a descriptive tumor-status module-trait comparison; the tumor-only analysis is the less status-dominated sensitivity network.",
  "- Macrophage-associated module correlations are reported in macrophage_module_association_all_samples.csv and macrophage_module_association_tumor_only.csv.",
  "",
  "## Highest target-set enrichment modules",
  "",
  "### All samples",
  "",
  if (nrow(all_top) == 0L) "No module enrichment rows were produced." else paste0("- ", all_top$module, ": ", all_top$target_gene_n_in_module, " target genes; Fisher P=", fmt(all_top$target_enrichment_p), "; BH-FDR=", fmt(all_top$target_enrichment_fdr), collapse = "\n"),
  "",
  "### Tumor-only",
  "",
  if (nrow(tumor_top) == 0L) "No module enrichment rows were produced." else paste0("- ", tumor_top$module, ": ", tumor_top$target_gene_n_in_module, " target genes; Fisher P=", fmt(tumor_top$target_enrichment_p), "; BH-FDR=", fmt(tumor_top$target_enrichment_fdr), collapse = "\n"),
  "",
  "## Interpretation boundary",
  "",
  "Target-set enrichment is calculated against the natural WGCNA input and uses only target genes selected by the prespecified MAD filter; all frozen targets remain visible in the overlay audit. Target/driver module membership or hubness should be treated as prioritization rather than mechanistic proof. Independent bulk or single-cell replication remains necessary."
)
writeLines(report, file.path(out_dir, "WGCNA_REPORT.md"), useBytes = TRUE)

manifest <- c(
  "analysis=dinp_crc_wgcna",
  "version=standardized_v2_unbiased_network_target_overlay",
  "dataset=GSE39582",
  "platform=GPL570",
  paste0("matrix_path=", normalizePath(matrix_path, winslash = "/", mustWork = FALSE)),
  paste0("annotation_path=", normalizePath(annotation_path, winslash = "/", mustWork = FALSE)),
  paste0("target_path=", normalizePath(target_path, winslash = "/", mustWork = FALSE)),
  paste0("driver_path=", normalizePath(driver_path, winslash = "/", mustWork = FALSE)),
  paste0("seed=", seed),
  paste0("max_genes=", max_genes),
  "targets_forced_into_wgcna_input=FALSE",
  paste0("macrophage_marker_set_n=", length(macrophage_markers)),
  "macrophage_score=within_subset_marker_z_mean_proxy",
  paste0("timestamp_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "network_type=signed",
  "correlation=bicor",
  "maxPOutliers=0.05",
  "minModuleSize=30",
  "mergeCutHeight=0.25",
  "probe_collapse=largest_across_sample_MAD",
  "all_sample_status_included=TRUE",
  "tumor_only_sensitivity=TRUE"
)
writeLines(manifest, file.path(out_dir, "WGCNA_MANIFEST.txt"), useBytes = TRUE)

message("WGCNA complete. Outputs written to: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
