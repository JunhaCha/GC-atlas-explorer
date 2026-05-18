#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(arrow)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  parsed <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected argument: %s", key))
    }
    key <- sub("^--", "", key)
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      parsed[[key]] <- TRUE
      i <- i + 1
    } else {
      parsed[[key]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  parsed
}

get_arg <- function(args, name, default = NULL, required = FALSE) {
  value <- args[[name]]
  if (is.null(value)) {
    if (required) stop(sprintf("Missing required argument --%s", name))
    return(default)
  }
  value
}

args <- parse_args()
input_path <- normalizePath(get_arg(args, "input", required = TRUE), winslash = "/", mustWork = TRUE)
bundle_path <- normalizePath(get_arg(args, "bundle", default = sub("\\.rds$", "_atlas_bundle_v2.rds", input_path, ignore.case = TRUE)), winslash = "/", mustWork = FALSE)
meta_path <- normalizePath(get_arg(args, "meta", default = sub("\\.rds$", "_atlas_meta_v2.parquet", input_path, ignore.case = TRUE)), winslash = "/", mustWork = FALSE)
umap_path <- normalizePath(get_arg(args, "umap", default = sub("\\.rds$", "_atlas_umap_v2.parquet", input_path, ignore.case = TRUE)), winslash = "/", mustWork = FALSE)
gene_dir <- normalizePath(get_arg(args, "gene-dir", default = sub("\\.rds$", "_atlas_gene_expr_v2", input_path, ignore.case = TRUE)), winslash = "/", mustWork = FALSE)
gene_manifest_path <- normalizePath(get_arg(args, "gene-manifest", default = file.path(gene_dir, "gene_manifest.rds")), winslash = "/", mustWork = FALSE)
assay_name <- get_arg(args, "assay", default = "RNA")
block_size <- as.integer(get_arg(args, "block-size", default = "128"))
if (!is.finite(block_size) || block_size <= 0) {
  stop("--block-size must be a positive integer.")
}

script_path <- normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]),
  winslash = "/",
  mustWork = FALSE
)
repo_dir <- dirname(script_path)
source(file.path(repo_dir, "R", "helpers.R"))

message("Reading Seurat object: ", input_path)
obj <- readRDS(input_path)
if (!inherits(obj, "Seurat")) {
  stop("Input is not a Seurat object.")
}

DefaultAssay(obj) <- assay_name
reduction_name <- preferred_umap_reduction(obj)
if (is.null(reduction_name)) {
  stop("No compatible UMAP reduction found in the Seurat object.")
}

coords <- Embeddings(obj, reduction = reduction_name)
coords <- as.data.frame(coords[, 1:2, drop = FALSE], stringsAsFactors = FALSE)
colnames(coords) <- c("UMAP_1", "UMAP_2")

meta <- obj[[]]
coords <- coords[rownames(meta), , drop = FALSE]
meta_out <- cbind(cell_id = rownames(meta), meta, stringsAsFactors = FALSE)
umap_out <- data.frame(
  cell_id = rownames(coords),
  UMAP_1 = coords$UMAP_1,
  UMAP_2 = coords$UMAP_2,
  stringsAsFactors = FALSE
)
expr <- get_expr_slot_safe(obj, assay_name = assay_name, slot_name = "data")
if (!inherits(expr, "dgCMatrix")) {
  expr <- as(expr, "dgCMatrix")
}

dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)

n_genes <- nrow(expr)
block_ids <- ((seq_len(n_genes) - 1L) %/% block_size) + 1L
gene_manifest <- data.frame(
  gene = rownames(expr),
  block_id = block_ids,
  row_in_block = ave(seq_len(n_genes), block_ids, FUN = seq_along),
  file_name = sprintf("block_%05d.rds", block_ids),
  stringsAsFactors = FALSE
)

message("Writing atlas expression block assets: ", gene_dir)
for (block_id in unique(block_ids)) {
  idx <- which(block_ids == block_id)
  block_mat <- expr[idx, , drop = FALSE]
  saveRDS(block_mat, file.path(gene_dir, sprintf("block_%05d.rds", block_id)), compress = TRUE)
}

message("Writing atlas gene manifest: ", gene_manifest_path)
saveRDS(gene_manifest, gene_manifest_path, compress = TRUE)

message("Writing atlas metadata parquet: ", meta_path)
arrow::write_parquet(meta_out, sink = meta_path)

message("Writing atlas UMAP parquet: ", umap_path)
arrow::write_parquet(umap_out, sink = umap_path)

bundle <- structure(
  list(
    type = "atlas_runtime_bundle",
    source_path = input_path,
    assay = assay_name,
    reduction_name = reduction_name,
    meta_path = meta_path,
    umap_path = umap_path,
    n_cells = nrow(meta),
    meta_colnames = colnames(meta),
    features = rownames(expr),
    gene_dir = gene_dir,
    gene_manifest_path = gene_manifest_path,
    gene_storage = "rds_sparse_block",
    gene_block_size = block_size
  ),
  class = "atlas_runtime_bundle"
)

message("Writing atlas runtime bundle: ", bundle_path)
saveRDS(bundle, bundle_path, compress = TRUE)
message("Done")
