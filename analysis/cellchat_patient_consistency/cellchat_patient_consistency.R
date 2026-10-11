#!/usr/bin/env Rscript

project_root <- "D:/DINP/whynot17-new-project"
out_dir <- file.path(project_root, "results", "2026-10-11", "cellchat_patient_consistency")
cellchat_lib <- file.path(project_root, "work", "Rlibrary")
.libPaths(unique(c(cellchat_lib, .libPaths())))
suppressPackageStartupMessages({
  library(CellChat)
  library(Matrix)
  library(dplyr)
  library(readr)
})

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = ""), "\n", sep = " ")
  flush.console()
}

read_cellchat_input <- function(cohort) {
  matrix_file <- file.path(out_dir, paste0(cohort, "_lognorm_genes_by_cells.mtx.gz"))
  genes_file <- file.path(out_dir, paste0(cohort, "_features.tsv"))
  cells_file <- file.path(out_dir, paste0(cohort, "_cells.tsv"))
  meta_file <- file.path(out_dir, paste0(cohort, "_cell_annotations.tsv.gz"))
  if (!file.exists(meta_file)) stop("Missing annotation file: ", meta_file)
  genes <- readLines(genes_file, warn = FALSE)
  cells <- readLines(cells_file, warn = FALSE)
  con <- gzfile(matrix_file, open = "rt")
  on.exit(close(con), add = TRUE)
  x <- Matrix::readMM(con)
  if (nrow(x) != length(genes) || ncol(x) != length(cells)) {
    stop(cohort, " matrix dimensions do not match features/cells sidecars")
  }
  rownames(x) <- make.unique(genes)
  colnames(x) <- cells
  x <- methods::as(x, "dgCMatrix")
  meta <- readr::read_tsv(meta_file, show_col_types = FALSE, progress = FALSE)
  if (!all(c("cell", "celltype", "sample_id", "cohort") %in% names(meta))) {
    stop(cohort, " metadata is missing required columns")
  }
  meta <- as.data.frame(meta)
  meta <- meta[match(cells, meta$cell), , drop = FALSE]
  if (anyNA(meta$cell)) stop(cohort, " metadata order does not match matrix cell names")
  rownames(meta) <- meta$cell
  list(data = x, meta = meta)
}

plus_one_bh <- function(tab, nboot) {
  if (is.null(tab) || nrow(tab) == 0L) return(tab)
  tab$pval_plus1 <- pmax(as.numeric(tab$pval), 1 / (nboot + 1))
  tab$qval_pair <- ave(tab$pval_plus1, interaction(tab$source, tab$target, drop = TRUE),
                       FUN = function(p) p.adjust(p, method = "BH"))
  tab
}

run_cohort <- function(cohort, min_cells = 10L, nboot = 100L) {
  inp <- read_cellchat_input(cohort)
  expr <- inp$data
  meta <- inp$meta
  meta$celltype <- as.character(meta$celltype)
  if (cohort == "GSE132465") {
    unit_col <- "patient_id"
    meta$unit_id <- as.character(meta$patient_id)
  } else {
    unit_col <- "sample_id"
    meta$unit_id <- as.character(meta$sample_id)
  }
  target_groups <- c("Epithelial", "Monocyte", "Macrophage", "Other_myeloid", "myCAF", "iCAF", "ECM_CAF", "apCAF")
  meta <- meta[meta$celltype %in% target_groups & !is.na(meta$unit_id) & nzchar(meta$unit_id), , drop = FALSE]
  expr <- expr[, meta$cell, drop = FALSE]
  caf_groups <- c("myCAF", "iCAF", "ECM_CAF", "apCAF")
  myel_groups <- c("Monocyte", "Macrophage", "Other_myeloid")
  units <- sort(unique(meta$unit_id))
  checkpoint_dir <- file.path(project_root, "work", "cellchat_checkpoints",
                              paste0(cohort, "_min", min_cells, "_nboot", nboot))
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  comm_rows <- list()
  elig_rows <- list()
  count_rows <- list()
  error_rows <- list()
  n_comm <- 0L
  n_elig <- 0L
  n_count <- 0L
  n_error <- 0L
  log_msg(cohort, ": starting ", length(units), if (cohort == "GSE132465") " patient units" else " tumor specimen units",
          "; CellChat v", as.character(packageVersion("CellChat")), "; nboot=", nboot)

  for (i in seq_along(units)) {
    unit <- units[[i]]
    md <- meta[meta$unit_id == unit, , drop = FALSE]
    group_n <- table(factor(md$celltype, levels = target_groups))
    group_n <- group_n[group_n > 0L]
    if (length(group_n) == 0L) next
    count_df <- data.frame(cohort = cohort, unit_id = unit, sample_id = unique(md$sample_id)[1],
                           celltype = names(group_n), n_cells = as.integer(group_n), stringsAsFactors = FALSE)
    n_count <- n_count + 1L
    count_rows[[n_count]] <- count_df
    # CellChat::filterCommunication(min.cells=10) removes groups with <=10 cells.
    # Match the eligibility denominator to the groups retained by that filter.
    keep_groups <- names(group_n)[group_n > min_cells]
    if (!("Epithelial" %in% keep_groups)) {
      log_msg(cohort, " / ", unit, ": excluded, fewer than ", min_cells, " epithelial cells")
      next
    }
    eligible_caf <- any(keep_groups %in% caf_groups)
    eligible_myel <- any(keep_groups %in% myel_groups)
    for (axis in c("Epithelial_CAF", "Epithelial_Myeloid")) {
      is_ok <- if (axis == "Epithelial_CAF") eligible_caf else eligible_myel
      if (is_ok) {
        n_elig <- n_elig + 1L
        elig_rows[[n_elig]] <- data.frame(cohort = cohort, unit_id = unit,
                                           sample_id = unique(md$sample_id)[1], axis = axis,
                                           epithelial_cells = as.integer(group_n[["Epithelial"]]),
                                           target_groups_present = paste(intersect(keep_groups, if (axis == "Epithelial_CAF") caf_groups else myel_groups), collapse = ";"),
                                           min_cells_per_group = min_cells,
                                           stringsAsFactors = FALSE)
      }
    }
    if (!(eligible_caf || eligible_myel)) {
      log_msg(cohort, " / ", unit, ": excluded, no CAF or myeloid group reached min.cells")
      next
    }
    checkpoint_file <- file.path(checkpoint_dir, paste0(gsub("[^A-Za-z0-9_-]", "_", unit), ".rds"))
    if (file.exists(checkpoint_file)) {
      cached_tab <- readRDS(checkpoint_file)
      cached_tab <- plus_one_bh(cached_tab, nboot)
      saveRDS(cached_tab, checkpoint_file)
      if (!is.null(cached_tab) && nrow(cached_tab) > 0L) {
        n_comm <- n_comm + 1L
        comm_rows[[n_comm]] <- cached_tab
      }
      log_msg(cohort, " / ", unit, ": reused patient/sample checkpoint")
      next
    }
    cell_keep <- md$cell[md$celltype %in% keep_groups]
    md_keep <- md[cell_keep, , drop = FALSE]
    x_keep <- expr[, cell_keep, drop = FALSE]
    md_keep$celltype <- factor(md_keep$celltype, levels = keep_groups)
    tryCatch({
      cc <- createCellChat(object = x_keep, meta = md_keep, group.by = "celltype", do.sparse = TRUE)
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
      if (!is.null(tab) && nrow(tab) > 0L) {
        tab$cohort <- cohort
        tab$unit_id <- unit
        tab$sample_id <- unique(md$sample_id)[1]
        tab$unit_type <- if (cohort == "GSE132465") "patient" else "tumor_sample"
        group_size <- data.frame(source = names(group_n), source_n = as.integer(group_n))
        tab <- tab %>% left_join(group_size, by = "source")
        group_size_target <- data.frame(target = names(group_n), target_n = as.integer(group_n))
        tab <- tab %>% left_join(group_size_target, by = "target")
        tab$axis <- dplyr::case_when(
          tab$source == "Epithelial" & tab$target %in% caf_groups ~ "Epithelial_CAF",
          tab$target == "Epithelial" & tab$source %in% caf_groups ~ "Epithelial_CAF",
          tab$source == "Epithelial" & tab$target %in% myel_groups ~ "Epithelial_Myeloid",
          tab$target == "Epithelial" & tab$source %in% myel_groups ~ "Epithelial_Myeloid",
          TRUE ~ NA_character_
        )
        tab$direction <- dplyr::case_when(
          tab$source == "Epithelial" ~ "Epithelial_to_target",
          tab$target == "Epithelial" ~ "Target_to_epithelial",
          TRUE ~ NA_character_
        )
        tab <- tab %>% filter(!is.na(axis))
        if (nrow(tab) > 0L) {
          tab$interaction_name <- if ("interaction_name" %in% names(tab)) as.character(tab$interaction_name) else ""
          tab$pathway_name <- if ("pathway_name" %in% names(tab)) as.character(tab$pathway_name) else ""
          tab$ligand_receptor <- paste(tab$ligand, tab$receptor, sep = "_")
          tab <- plus_one_bh(tab, nboot)
          n_comm <- n_comm + 1L
          comm_rows[[n_comm]] <- tab
        }
      }
      saveRDS(tab, checkpoint_file)
      rm(cc)
      log_msg(cohort, " / ", unit, ": completed; groups=", paste(names(group_n), as.integer(group_n), sep = ":", collapse = ","),
              "; retained for axes: CAF=", eligible_caf, ", myeloid=", eligible_myel)
    }, error = function(e) {
      n_error <<- n_error + 1L
      error_rows[[n_error]] <<- data.frame(cohort = cohort, unit_id = unit, sample_id = unique(md$sample_id)[1],
                                           error = conditionMessage(e), stringsAsFactors = FALSE)
      log_msg(cohort, " / ", unit, ": ERROR ", conditionMessage(e))
    })
    rm(x_keep, md_keep)
    gc(verbose = FALSE)
  }
  comm <- if (length(comm_rows)) bind_rows(comm_rows) else data.frame()
  eligible <- if (length(elig_rows)) bind_rows(elig_rows) else data.frame()
  cell_counts <- if (length(count_rows)) bind_rows(count_rows) else data.frame()
  errors <- if (length(error_rows)) bind_rows(error_rows) else data.frame()
  if (nrow(cell_counts)) readr::write_tsv(cell_counts, file.path(out_dir, paste0(cohort, "_CellChat_group_cell_counts.tsv")))
  if (nrow(eligible)) readr::write_tsv(eligible, file.path(out_dir, paste0(cohort, "_CellChat_axis_eligible_units.tsv")))
  if (nrow(errors)) readr::write_tsv(errors, file.path(out_dir, paste0(cohort, "_CellChat_errors.tsv")))
  if (!nrow(errors)) writeLines("cohort\tunit_id\tsample_id\terror", file.path(out_dir, paste0(cohort, "_CellChat_errors.tsv")))
  if (nrow(comm)) {
    readr::write_tsv(comm, file.path(out_dir, paste0(cohort, "_CellChat_patient_sample_interactions.tsv.gz")))
  }
  list(comm = comm, eligible = eligible, counts = cell_counts, errors = errors)
}

summarize_cohort <- function(x, cohort) {
  tab <- x$comm
  eligible <- x$eligible
  if (is.null(tab) || !nrow(tab)) return(data.frame())
  edge <- tab %>%
    filter(prob > 0) %>%
    group_by(cohort, axis, direction, source, target, ligand, receptor, ligand_receptor, pathway_name, interaction_name) %>%
    summarise(
      n_detected = n_distinct(unit_id),
      n_raw_p05 = n_distinct(unit_id[pval_plus1 < 0.05]),
      n_pair_q05 = n_distinct(unit_id[qval_pair < 0.05]),
      median_prob = median(prob, na.rm = TRUE),
      mean_prob = mean(prob, na.rm = TRUE),
      median_raw_p = median(pval_plus1, na.rm = TRUE),
      median_pair_q = median(qval_pair, na.rm = TRUE),
      .groups = "drop"
    )
  denominators <- eligible %>% group_by(cohort, axis) %>% summarise(n_eligible_axis = n_distinct(unit_id), .groups = "drop")
  group_counts <- x$counts %>% filter(n_cells > 10L)
  pair_denominators <- inner_join(group_counts, group_counts,
                                  by = c("cohort", "unit_id", "sample_id"), suffix = c("_source", "_target")) %>%
    transmute(cohort, unit_id, source = celltype_source, target = celltype_target,
              axis = case_when(
                source == "Epithelial" & target %in% c("myCAF", "iCAF", "ECM_CAF", "apCAF") ~ "Epithelial_CAF",
                target == "Epithelial" & source %in% c("myCAF", "iCAF", "ECM_CAF", "apCAF") ~ "Epithelial_CAF",
                source == "Epithelial" & target %in% c("Monocyte", "Macrophage", "Other_myeloid") ~ "Epithelial_Myeloid",
                target == "Epithelial" & source %in% c("Monocyte", "Macrophage", "Other_myeloid") ~ "Epithelial_Myeloid",
                TRUE ~ NA_character_
              )) %>% filter(!is.na(axis)) %>%
    group_by(cohort, axis, source, target) %>% summarise(n_eligible_pair = n_distinct(unit_id), .groups = "drop")
  edge %>%
    left_join(denominators, by = c("cohort", "axis")) %>%
    left_join(pair_denominators, by = c("cohort", "axis", "source", "target")) %>%
    mutate(raw_p05_fraction = n_raw_p05 / pmax(n_eligible_pair, 1),
           pair_q05_fraction = n_pair_q05 / pmax(n_eligible_pair, 1),
           detection_fraction = n_detected / pmax(n_eligible_pair, 1),
           unit_type = ifelse(cohort == "GSE132465", "patient", "tumor_sample"))
}

main <- function() {
  log_msg("CellChat version ", as.character(packageVersion("CellChat")),
          "; human database interactions=", nrow(CellChatDB.human$interaction),
          "; min.cells=10; nboot=100; population.size=FALSE")
  res132 <- run_cohort("GSE132465", min_cells = 10L, nboot = 100L)
  s132 <- summarize_cohort(res132, "GSE132465")
  if (nrow(s132)) readr::write_tsv(s132, file.path(out_dir, "GSE132465_CellChat_recurrence_by_patient.tsv.gz"))
  res132$comm <- NULL
  gc(verbose = FALSE)
  res200 <- run_cohort("GSE200997", min_cells = 10L, nboot = 100L)
  s200 <- summarize_cohort(res200, "GSE200997")
  if (nrow(s200)) readr::write_tsv(s200, file.path(out_dir, "GSE200997_CellChat_recurrence_by_sample.tsv.gz"))
  if (nrow(s132) && nrow(s200)) {
    keys <- c("axis", "direction", "source", "target", "ligand", "receptor", "ligand_receptor", "pathway_name", "interaction_name")
    combined <- full_join(
      s132 %>% select(all_of(keys), GSE132465_n_eligible = n_eligible_pair,
                      GSE132465_n_patient_raw_p05 = n_raw_p05, GSE132465_patient_raw_p05_fraction = raw_p05_fraction,
                      GSE132465_n_patient_q05 = n_pair_q05, GSE132465_patient_q05_fraction = pair_q05_fraction,
                      GSE132465_median_prob = median_prob),
      s200 %>% select(all_of(keys), GSE200997_n_eligible = n_eligible_pair,
                      GSE200997_n_sample_raw_p05 = n_raw_p05, GSE200997_sample_raw_p05_fraction = raw_p05_fraction,
                      GSE200997_n_sample_q05 = n_pair_q05, GSE200997_sample_q05_fraction = pair_q05_fraction,
                      GSE200997_median_prob = median_prob), by = keys
    ) %>% mutate(
      across(starts_with("GSE"), ~ tidyr::replace_na(.x, 0)),
      replicated_raw_p05_at_least_half = GSE132465_n_eligible >= 2 & GSE200997_n_eligible >= 2 &
        GSE132465_patient_raw_p05_fraction >= 0.5 & GSE200997_sample_raw_p05_fraction >= 0.5,
      replicated_pair_q05_at_least_half = GSE132465_n_eligible >= 2 & GSE200997_n_eligible >= 2 &
        GSE132465_patient_q05_fraction >= 0.5 & GSE200997_sample_q05_fraction >= 0.5
    ) %>% arrange(desc(replicated_pair_q05_at_least_half), desc(replicated_raw_p05_at_least_half),
                  desc(GSE132465_patient_raw_p05_fraction), desc(GSE200997_sample_raw_p05_fraction), desc(GSE132465_median_prob))
    readr::write_tsv(combined, file.path(out_dir, "CellChat_cross_cohort_sensitivity.tsv.gz"))
  }
  eligibility <- bind_rows(res132$eligible, res200$eligible)
  readr::write_tsv(eligibility, file.path(out_dir, "CellChat_axis_eligible_units_all_cohorts.tsv"))
  summary <- eligibility %>% group_by(cohort, axis) %>% summarise(n_eligible_units = n_distinct(unit_id), .groups = "drop")
  readr::write_tsv(summary, file.path(out_dir, "CellChat_axis_eligible_unit_summary.tsv"))
  log_msg("CellChat patient/sample analysis finished")
  print(summary)
}

main()
