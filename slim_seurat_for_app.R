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

split_csv <- function(x) {
  parts <- trimws(unlist(strsplit(x, ",")))
  parts[nzchar(parts)]
}

args <- parse_args()

input_path <- get_arg(args, "input", required = TRUE)
output_path <- get_arg(args, "output", required = TRUE)
assay_name <- get_arg(args, "assay", default = NULL)
reduction_names <- split_csv(get_arg(args, "reductions", default = "HarmonyUMAP"))
layer_names <- split_csv(get_arg(args, "layers", default = "data"))
keep_misc <- isTRUE(as.logical(get_arg(args, "keep_misc", default = "FALSE")))

if (!file.exists(input_path)) {
  stop("Input Seurat object not found: ", input_path)
}

message("Reading: ", input_path)
obj <- readRDS(input_path)
if (!inherits(obj, "Seurat")) {
  stop("Input file is not a Seurat object.")
}

assay_to_use <- assay_name %||% DefaultAssay(obj)
if (!assay_to_use %in% names(obj@assays)) {
  stop("Assay not found in object: ", assay_to_use)
}

reductions_to_keep <- reduction_names[reduction_names %in% names(obj@reductions)]
if (length(reductions_to_keep) == 0) {
  warning("None of the requested reductions were found. The slim object will keep no reductions.")
}

message("Keeping assay: ", assay_to_use)
message("Keeping layers: ", paste(layer_names, collapse = ", "))
message("Keeping reductions: ", if (length(reductions_to_keep) > 0) paste(reductions_to_keep, collapse = ", ") else "<none>")

slim <- DietSeurat(
  object = obj,
  assays = assay_to_use,
  layers = layer_names,
  dimreducs = reductions_to_keep,
  graphs = NULL,
  misc = keep_misc
)

DefaultAssay(slim) <- assay_to_use

if ("commands" %in% slotNames(slim)) {
  slot(slim, "commands") <- list()
}
if ("tools" %in% slotNames(slim)) {
  slot(slim, "tools") <- list()
}

gc()

message("Writing: ", output_path)
saveRDS(slim, output_path, compress = TRUE)

input_size <- file.info(input_path)$size / 1024^3
output_size <- file.info(output_path)$size / 1024^3
message(sprintf("Input size: %.2f GB", input_size))
message(sprintf("Output size: %.2f GB", output_size))
