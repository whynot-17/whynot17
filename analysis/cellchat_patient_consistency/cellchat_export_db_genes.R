project_root <- "D:/DINP/whynot17-new-project"
lib_cellchat <- file.path(project_root, "work", "Rlibrary")
.libPaths(unique(c(lib_cellchat, .libPaths())))
suppressPackageStartupMessages(library(CellChat))

db <- CellChatDB.human$interaction
raw_symbols <- unique(c(db$ligand.symbol, db$receptor.symbol,
                        db$ligand, db$receptor))
genes <- unique(unlist(strsplit(raw_symbols, "[_+|,; ]+", perl = TRUE)))
genes <- sort(genes[nzchar(genes) & !is.na(genes)])
out <- file.path(project_root, "work", "cellchat_db_gene_symbols.txt")
writeLines(genes, out, useBytes = TRUE)
cat("CellChat version:", as.character(packageVersion("CellChat")), "\n")
cat("Interaction database rows:", nrow(db), "\n")
cat("Unique ligand/receptor symbols:", length(genes), "\n")
cat("Wrote:", out, "\n")
