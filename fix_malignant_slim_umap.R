#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) {
    y
  } else {
    x
  }
}

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
    if (required) {
      stop(sprintf("Missing required argument --%s", name))
    }
    return(default)
  }
  value
}

args <- parse_args()
raw_path <- normalizePath(get_arg(args, "raw", required = TRUE), winslash = "/", mustWork = TRUE)
slim_path <- normalizePath(get_arg(args, "slim", required = TRUE), winslash = "/", mustWork = TRUE)
reduction_name <- get_arg(args, "reduction", default = "umap")

message("Reading raw object: ", raw_path)
raw_obj <- readRDS(raw_path)
message("Reading slim object: ", slim_path)
slim_obj <- readRDS(slim_path)

if (!inherits(raw_obj, "Seurat") || !inherits(slim_obj, "Seurat")) {
  stop("Both input files must be Seurat objects.")
}

if (!reduction_name %in% names(raw_obj@reductions)) {
  stop("Reduction not found in raw object: ", reduction_name)
}

if (!identical(colnames(raw_obj), colnames(slim_obj))) {
  stop("Raw and slim object cells are not in the same order; refusing to patch.")
}

message("Copying reduction '", reduction_name, "' into slim object.")
slim_obj[[reduction_name]] <- raw_obj[[reduction_name]]

message("Writing patched slim object: ", slim_path)
saveRDS(slim_obj, slim_path, compress = TRUE)
message("Done")
