#!/usr/bin/env Rscript

# Re-run patient/sample-level CellChat for tumor epithelium and the current
# PID1-independent myeloid subtype labels. The old normalized ligand/receptor
# panel matrix and epithelial-only labels are reused; old myeloid labels and
# old communication results are not used as model inputs.

project_root <- "D:/DINP/whynot17-new-project"
old_dir <- file.path(project_root, "results", "2026-10-11", "cellchat_patient_consistency")
ann_dir <- file.path(project_root, "results", "2026-10-11", "PID1_myeloid_finelocalization")
out_dir <- file.path(project_root, "results", "2026-10-11", "PID1_CellChat_subtype_consistency")
checkpoint_root <- file.path(project_root, "work", "cellchat_pid1_subtype_checkpoints_20261011")
cellchat_lib <- file.path(project_root, "work", "Rlibrary")
.libPaths(unique(c(cellchat_lib, .libPaths())))
suppressPackageStartupMessages({
  library(CellChat)
  library(Matrix)
  library(dplyr)
  library(readr)
  library(tidyr)
})

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_root, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = ""), "\n", sep = " ")
  flush.console()
}

read_input <- function(cohort) {
  matrix_file <- file.path(old_dir, paste0(cohort, "_lognorm_genes_by_cells.mtx.gz"))
  genes_file <- file.path(old_dir, paste0(cohort, "_features.tsv"))
  cells_file <- file.path(old_dir, paste0(cohort, "_cells.tsv"))
  old_meta_file <- file.path(old_dir, paste0(cohort, "_cell_annotations.tsv.gz"))
  new_meta_file <- file.path(ann_dir, paste0(cohort, "_PID1_CellChat_myeloid_annotations.tsv.gz"))

  genes <- readLines(genes_file, warn = FALSE)
  cells <- readLines(cells_file, warn = FALSE)
  con <- gzfile(matrix_file, open = "rt")
  on.exit(close(con), add = TRUE)
  expr <- Matrix::readMM(con)
  if (nrow(expr) != length(genes) || ncol(expr) != length(cells)) {
    stop(cohort, ": expression matrix dimensions disagree with sidecars")
  }
  rownames(expr) <- make.unique(genes)
  colnames(expr) <- cells
  expr <- methods::as(expr, "dgCMatrix")

  old_meta <- readr::read_tsv(old_meta_file, show_col_types = FALSE, progress = FALSE) %>% as.data.frame()
  cur_meta <- readr::read_tsv(new_meta_file, show_col_types = FALSE, progress = FALSE) %>% as.data.frame()
  if (!all(c("cell", "celltype", "sample_id", "patient_id") %in% names(old_meta))) {
    stop(cohort, ": old epithelial annotation lacks required fields")
  }
  if (!all(c("cell_id", "unit_id", "pid1_subtype", "cellchat_eligible", "PID1_log1p", "PID1_detected") %in% names(cur_meta))) {
    stop(cohort, ": current PID1 annotation lacks required fields")
  }
  if (cohort == "GSE132465") {
    epi_unit <- as.character(old_meta$patient_id)
    unit_type <- "patient"
  } else {
    epi_unit <- as.character(old_meta$sample_id)
    unit_type <- "tumor_sample"
  }
  epi_keep <- old_meta$celltype == "Epithelial" & !is.na(epi_unit) & nzchar(epi_unit)
  epi <- data.frame(
    cell = as.character(old_meta$cell[epi_keep]),
    sample_id = as.character(old_meta$sample_id[epi_keep]),
    unit_id = epi_unit[epi_keep],
    celltype = "Epithelial",
    pid1_subtype = NA_character_,
    PID1_log1p = NA_real_,
    PID1_detected = NA,
    cohort = cohort,
    unit_type = unit_type,
    stringsAsFactors = FALSE
  )
  is_eligible <- tolower(as.character(cur_meta$cellchat_eligible)) %in% c("true", "1", "yes")
  is_mono_mac <- grepl("^(Monocyte|Macrophage)", as.character(cur_meta$pid1_subtype))
  my <- cur_meta[is_eligible & is_mono_mac, , drop = FALSE]
  if (cohort == "GSE132465") {
    my_unit <- as.character(my$unit_id)
  } else {
    my_unit <- as.character(my$sample_id)
  }
  my <- data.frame(
    cell = as.character(my$cell_id),
    sample_id = as.character(my$sample_id),
    unit_id = my_unit,
    celltype = as.character(my$pid1_subtype),
    pid1_subtype = as.character(my$pid1_subtype),
    PID1_log1p = as.numeric(my$PID1_log1p),
    PID1_detected = as.logical(my$PID1_detected),
    cohort = cohort,
    unit_type = unit_type,
    stringsAsFactors = FALSE
  )
  meta <- bind_rows(epi, my)
  if (anyDuplicated(meta$cell)) stop(cohort, ": duplicate cells after epithelial/myeloid merge")
  missing_cells <- setdiff(meta$cell, colnames(expr))
  if (length(missing_cells)) stop(cohort, ": current cells missing from CellChat expression matrix: ", length(missing_cells))
  meta <- meta[match(colnames(expr)[colnames(expr) %in% meta$cell], meta$cell), , drop = FALSE]
  expr <- expr[, meta$cell, drop = FALSE]
  if (!identical(colnames(expr), meta$cell)) stop(cohort, ": expression and metadata order mismatch")
  rownames(meta) <- meta$cell

  # Unit/subtype expression and count table used for eligibility and PID1 association.
  cell_counts <- meta %>%
    group_by(cohort, unit_id, sample_id, unit_type, celltype, pid1_subtype, .drop = TRUE) %>%
    summarise(
      n_cells = dplyr::n(),
      PID1_mean_log1p = if (all(is.na(PID1_log1p))) NA_real_ else mean(PID1_log1p, na.rm = TRUE),
      PID1_detection_fraction = if (all(is.na(PID1_detected))) NA_real_ else mean(PID1_detected, na.rm = TRUE),
      .groups = "drop"
    )
  list(expr = expr, meta = meta, cell_counts = cell_counts)
}

plus_one_bh <- function(tab, nboot) {
  if (is.null(tab) || nrow(tab) == 0L) return(tab)
  p <- suppressWarnings(as.numeric(tab$pval))
  p[is.na(p)] <- 1
  tab$pval_plus1 <- pmax(p, 1 / (nboot + 1))
  tab$qval_pair <- ave(tab$pval_plus1, interaction(tab$source, tab$target, drop = TRUE),
                       FUN = function(x) p.adjust(x, method = "BH"))
  tab
}

run_cohort <- function(cohort, min_cells = 10L, nboot = 100L) {
  inp <- read_input(cohort)
  expr <- inp$expr
  meta <- inp$meta
  counts <- inp$cell_counts
  groups <- c(
    "Epithelial", "Monocyte_CD14_classical", "Monocyte_CD14_classical_inflammatory",
    "Macrophage_C1QC", "Macrophage_inflammatory", "Macrophage_SPP1_TREM2", "Macrophage_IFN"
  )
  units <- sort(unique(meta$unit_id))
  pair_rows <- list()
  comm_rows <- list()
  count_index <- 0L
  comm_index <- 0L
  elig_index <- 0L
  elig_rows <- list()
  checkpoint_dir <- file.path(checkpoint_root, paste0(cohort, "_min", min_cells, "_nboot", nboot))
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

  log_msg(cohort, ": cells=", ncol(expr), "; units=", length(units), "; min.cells=", min_cells,
          "; nboot=", nboot, "; CellChat=", as.character(packageVersion("CellChat")))
  for (i in seq_along(units)) {
    unit <- units[[i]]
    md <- meta[meta$unit_id == unit, , drop = FALSE]
    n_by <- table(factor(md$celltype, levels = groups))
    n_by <- n_by[n_by > 0L]
    n_epi <- as.integer(n_by[["Epithelial"]] %||% 0L)
    if (is.na(n_epi) || n_epi <= min_cells) next
    target_groups <- setdiff(names(n_by), "Epithelial")
    target_groups <- target_groups[as.integer(n_by[target_groups]) > min_cells]
    if (!length(target_groups)) next
    for (state in target_groups) {
      cc_row <- counts %>% filter(unit_id == unit, celltype == state)
      elig_index <- elig_index + 1L
      elig_rows[[elig_index]] <- data.frame(
        cohort = cohort, unit_id = unit, sample_id = unique(md$sample_id)[1], unit_type = unique(md$unit_type)[1],
        pid1_subtype = state, n_epithelial_cells = n_epi, n_myeloid_cells = as.integer(n_by[[state]]),
        PID1_mean_log1p = cc_row$PID1_mean_log1p[1],
        PID1_detection_fraction = cc_row$PID1_detection_fraction[1],
        eligible_gt10 = as.integer(n_by[[state]]) > 10L,
        eligible_ge20 = as.integer(n_by[[state]]) >= 20L,
        eligible_ge50 = as.integer(n_by[[state]]) >= 50L,
        stringsAsFactors = FALSE
      )
    }
    keep_groups <- c("Epithelial", target_groups)
    cell_keep <- md$cell[md$celltype %in% keep_groups]
    md_keep <- md[match(cell_keep, md$cell), , drop = FALSE]
    md_keep$celltype <- factor(md_keep$celltype, levels = keep_groups)
    expr_keep <- expr[, cell_keep, drop = FALSE]
    checkpoint_file <- file.path(checkpoint_dir, paste0(gsub("[^A-Za-z0-9_-]", "_", unit), ".rds"))
    if (file.exists(checkpoint_file)) {
      tab <- readRDS(checkpoint_file)
      if (is.null(tab) || !(".error" %in% names(tab))) {
        if (!is.null(tab) && nrow(tab)) {
          comm_index <- comm_index + 1L
          comm_rows[[comm_index]] <- tab
        }
        log_msg(cohort, " / ", unit, ": reused this-run checkpoint")
        next
      }
      log_msg(cohort, " / ", unit, ": retrying a previously failed this-run checkpoint")
      unlink(checkpoint_file)
    }
    tryCatch({
      cc <- createCellChat(object = expr_keep, meta = md_keep, group.by = "celltype", do.sparse = TRUE)
      cc@DB <- CellChatDB.human
      cc <- subsetData(cc)
      cc <- identifyOverExpressedGenes(cc, thresh.p = 0.05, min.cells = min_cells, do.fast = FALSE)
      cc <- identifyOverExpressedInteractions(cc)
      cc <- computeCommunProb(cc, type = "triMean", raw.use = TRUE, population.size = FALSE,
                              nboot = nboot, seed.use = 20261011 + i)
      cc <- filterCommunication(cc, min.cells = min_cells)
      cc <- computeCommunProbPathway(cc)
      cc <- aggregateNet(cc, thresh = 1)
      tab <- subsetCommunication(cc, slot.name = "net", thresh = 1)
      if (is.null(tab)) tab <- data.frame()
      if (nrow(tab)) {
        tab$source <- as.character(tab$source)
        tab$target <- as.character(tab$target)
        tab$cohort <- cohort
        tab$unit_id <- unit
        tab$sample_id <- unique(md$sample_id)[1]
        tab$unit_type <- unique(md$unit_type)[1]
        tab$direction <- dplyr::case_when(
          tab$source == "Epithelial" & tab$target != "Epithelial" ~ "Epithelial_to_myeloid",
          tab$target == "Epithelial" & tab$source != "Epithelial" ~ "Myeloid_to_epithelial",
          TRUE ~ NA_character_
        )
        tab$pid1_subtype <- as.character(ifelse(tab$direction == "Epithelial_to_myeloid", tab$target, tab$source))
        tab <- tab %>% filter(!is.na(direction))
        tab$ligand_receptor <- paste(as.character(tab$ligand), as.character(tab$receptor), sep = "_")
        tab$interaction_name <- if ("interaction_name" %in% names(tab)) as.character(tab$interaction_name) else tab$ligand_receptor
        # subsetCommunication can return rows that are not epithelial–myeloid
        # in a sparse sample; after direction filtering the table may be empty.
        # Use a zero-length vector then so a no-edge patient is recorded as a
        # successful zero-edge model rather than an R scalar-assignment error.
        tab$n_epithelial_cells <- rep(n_epi, nrow(tab))
        tab$n_myeloid_cells <- as.integer(n_by[tab$pid1_subtype])
        tab <- tab %>% left_join(
          counts %>% filter(celltype != "Epithelial") %>%
            select(cohort, unit_id, pid1_subtype = celltype, PID1_mean_log1p, PID1_detection_fraction),
          by = c("cohort", "unit_id", "pid1_subtype")
        )
        tab <- plus_one_bh(tab, nboot)
        comm_index <- comm_index + 1L
        comm_rows[[comm_index]] <- tab
      }
      saveRDS(tab, checkpoint_file)
      rm(cc)
      log_msg(cohort, " / ", unit, ": completed; groups=", paste(names(n_by), as.integer(n_by), sep = ":", collapse = ","))
    }, error = function(e) {
      log_msg(cohort, " / ", unit, ": ERROR ", conditionMessage(e))
      saveRDS(data.frame(.error = conditionMessage(e)), checkpoint_file)
    })
    rm(expr_keep, md_keep)
    gc(verbose = FALSE)
  }
  comm <- if (length(comm_rows)) bind_rows(comm_rows) else data.frame()
  elig <- if (length(elig_rows)) bind_rows(elig_rows) else data.frame()
  # A checkpoint can be a successful zero-row table or an error record; expose any errors.
  checkpoint_files <- list.files(checkpoint_dir, pattern = "\\.rds$", full.names = TRUE)
  errors <- lapply(checkpoint_files, function(p) {
    z <- tryCatch(readRDS(p), error = function(e) NULL)
    if (!is.null(z) && ".error" %in% names(z)) data.frame(cohort = cohort, unit_id = basename(p), error = z$.error[1]) else NULL
  }) %>% bind_rows()
  if (!nrow(errors)) errors <- data.frame(cohort = character(), unit_id = character(), error = character())
  readr::write_tsv(inp$cell_counts, file.path(out_dir, paste0(cohort, "_PID1_CellChat_cell_counts_all_groups.tsv")))
  readr::write_tsv(elig, file.path(out_dir, paste0(cohort, "_PID1_CellChat_pair_eligibility.tsv")))
  readr::write_tsv(errors, file.path(out_dir, paste0(cohort, "_PID1_CellChat_errors.tsv")))
  readr::write_tsv(meta, file.path(out_dir, paste0(cohort, "_PID1_CellChat_epithelial_myeloid_cell_annotations.tsv.gz")))
  if (nrow(comm)) readr::write_tsv(comm, file.path(out_dir, paste0(cohort, "_PID1_CellChat_patient_sample_interactions.tsv.gz")))
  list(comm = comm, eligible = elig, counts = inp$cell_counts, errors = errors)
}

main <- function() {
  log_msg("Starting PID1 subtype CellChat analysis; reusing only previous CellChat panel matrices and epithelial labels, with current myeloid subtype annotations")
  log_msg("CellChat=", as.character(packageVersion("CellChat")), "; DB interactions=", nrow(CellChatDB.human$interaction),
          "; min.cells=10; nboot=100; population.size=FALSE")
  a <- run_cohort("GSE132465", min_cells = 10L, nboot = 100L)
  gc(verbose = FALSE)
  b <- run_cohort("GSE200997", min_cells = 10L, nboot = 100L)
  all_elig <- bind_rows(a$eligible, b$eligible)
  all_comm <- bind_rows(a$comm, b$comm)
  readr::write_tsv(all_elig, file.path(out_dir, "PID1_CellChat_pair_eligibility_all_cohorts.tsv"))
  if (nrow(all_comm)) readr::write_tsv(all_comm, file.path(out_dir, "PID1_CellChat_patient_sample_interactions_all_cohorts.tsv.gz"))
  summary <- all_elig %>% group_by(cohort, pid1_subtype) %>% summarise(
    n_eligible_gt10 = n_distinct(unit_id),
    n_eligible_ge20 = n_distinct(unit_id[eligible_ge20]),
    n_eligible_ge50 = n_distinct(unit_id[eligible_ge50]),
    median_subtype_cells = median(n_myeloid_cells),
    median_PID1_mean_log1p = median(PID1_mean_log1p, na.rm = TRUE),
    median_PID1_detection = median(PID1_detection_fraction, na.rm = TRUE),
    .groups = "drop"
  )
  readr::write_tsv(summary, file.path(out_dir, "PID1_CellChat_subtype_coverage_summary.tsv"))
  log_msg("Finished; eligible pairs=", nrow(all_elig), "; interaction rows=", nrow(all_comm),
          "; errors=", nrow(a$errors) + nrow(b$errors))
  print(summary)
}

`%||%` <- function(x, y) if (length(x) == 0L || is.na(x)) y else x
main()
