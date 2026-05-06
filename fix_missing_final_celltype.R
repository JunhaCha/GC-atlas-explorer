#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1 && nzchar(args[[1]])) {
  args[[1]]
} else {
  "/Users/junhacha/Documents/Playground"
}

specific_file <- if (length(args) >= 2 && nzchar(args[[2]])) args[[2]] else NULL

if (!is.null(specific_file)) {
  files <- normalizePath(specific_file, winslash = "/", mustWork = FALSE)
} else {
  files <- sort(Sys.glob(file.path(data_dir, "seurat*.rds")))
  files <- files[!(endsWith(files, "_avg_cache.rds") | endsWith(files, "_quadrant_cache.rds"))]
}

updated <- character(0)
skipped <- character(0)

for (f in files) {
  cat("[check]", basename(f), "\n")
  obj <- readRDS(f)

  if (!inherits(obj, "Seurat")) {
    cat("  not a Seurat object, skip\n")
    skipped <- c(skipped, basename(f))
    next
  }

  md <- obj@meta.data
  if ("final_celltype" %in% colnames(md)) {
    cat("  final_celltype already present\n")
    skipped <- c(skipped, basename(f))
    next
  }

  md$final_celltype <- as.character(Idents(obj))
  obj@meta.data <- md
  saveRDS(obj, f, compress = TRUE)

  cat("  added final_celltype from active.ident and saved\n")
  updated <- c(updated, basename(f))
}

cat("\nUPDATED\n")
if (length(updated)) {
  cat(paste(updated, collapse = "\n"), "\n")
} else {
  cat("none\n")
}

cat("\nSKIPPED\n")
if (length(skipped)) {
  cat(paste(skipped, collapse = "\n"), "\n")
} else {
  cat("none\n")
}
