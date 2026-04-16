#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
})

source("R/helpers.R")

safe_file_stub <- function(x) {
  gsub("[^A-Za-z0-9._-]", "_", x)
}

build_sample_marker_assets <- function(sample_id) {
  counts_path <- xenium_counts_matrix_path(sample_id)
  cells_path <- file.path(xenium_app_inputs_dir(), sample_id, "cells.rds")

  if (!file.exists(counts_path)) {
    stop(paste0("Counts matrix not found for sample ", sample_id, ": ", counts_path), call. = FALSE)
  }
  if (!file.exists(cells_path)) {
    stop(paste0("Cell table not found for sample ", sample_id, ": ", cells_path), call. = FALSE)
  }

  message("[xenium-markers] Loading sample ", sample_id)
  mat <- readRDS(counts_path)
  if (!inherits(mat, "dgCMatrix")) {
    stop(paste0("Expected a dgCMatrix for sample ", sample_id, "."), call. = FALSE)
  }

  genes_all <- mat@Dimnames[[1]]
  cells_all <- mat@Dimnames[[2]]
  cells_df <- readRDS(cells_path)

  if (is.null(genes_all) || is.null(cells_all)) {
    stop(paste0("Counts matrix dimnames are missing for sample ", sample_id, "."), call. = FALSE)
  }

  cell_match <- match(cells_df$cell_id, cells_all)
  if (anyNA(cell_match)) {
    stop(paste0("Cell IDs in cells.rds do not fully match the counts matrix for sample ", sample_id, "."), call. = FALSE)
  }
  if (!identical(cell_match, seq_along(cells_all))) {
    message("[xenium-markers] Reordering matrix columns to match cells.rds for ", sample_id)
    mat <- mat[, cell_match, drop = FALSE]
  }

  genes_use <- genes_all[!is.na(genes_all) & nzchar(genes_all)]
  gene_idx <- match(genes_use, genes_all)
  sub_mat <- mat[gene_idx, , drop = FALSE]
  trip <- Matrix::summary(sub_mat)

  out_dir <- xenium_marker_dir(sample_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  old_rds <- list.files(out_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(old_rds) > 0) {
    unlink(old_rds)
  }

  split_rows <- split(seq_len(nrow(trip)), trip$i)
  manifest_rows <- vector("list", length(genes_use))

  for (i in seq_along(genes_use)) {
    gene <- genes_use[[i]]
    file_name <- paste0(safe_file_stub(gene), ".rds")
    rows <- split_rows[[as.character(i)]]

    if (length(rows) > 0) {
      marker_obj <- list(
        gene = gene,
        index = as.integer(trip$j[rows]),
        value = as.numeric(trip$x[rows]),
        n_cells = ncol(sub_mat)
      )
    } else {
      marker_obj <- list(
        gene = gene,
        index = integer(0),
        value = numeric(0),
        n_cells = ncol(sub_mat)
      )
    }

    saveRDS(marker_obj, file.path(out_dir, file_name))
    manifest_rows[[i]] <- data.frame(
      gene = gene,
      file_name = file_name,
      stringsAsFactors = FALSE
    )
  }

  manifest <- do.call(rbind, manifest_rows)
  saveRDS(manifest, file.path(out_dir, "marker_genes.rds"))

  message(
    "[xenium-markers] ",
    sample_id,
    ": wrote ",
    nrow(manifest),
    " marker files to ",
    out_dir
  )
}

main <- function() {
  manifest <- xenium_sample_manifest()

  for (sample_id in manifest$sample_id) {
    build_sample_marker_assets(sample_id)
  }

  message("[xenium-markers] Finished building marker assets for all samples.")
}

main()
