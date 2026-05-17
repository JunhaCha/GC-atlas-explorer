#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
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
    if (required) {
      stop(sprintf("Missing required argument --%s", name))
    }
    return(default)
  }
  value
}

args <- parse_args()
input_path <- normalizePath(get_arg(args, "input", required = TRUE), winslash = "/", mustWork = TRUE)
output_path <- normalizePath(get_arg(args, "output", default = input_path), winslash = "/", mustWork = FALSE)
lineage_label <- get_arg(args, "label", default = "malignant")

message("Reading: ", input_path)
obj <- readRDS(input_path)
if (!inherits(obj, "Seurat")) {
  stop("Input file is not a Seurat object.")
}

md <- obj@meta.data
if (!"rev_pathological_subtype" %in% colnames(md)) {
  stop("rev_pathological_subtype column is required.")
}

keep_cells <- rownames(md)[!is.na(md$rev_pathological_subtype) & as.character(md$rev_pathological_subtype) != "Normal"]
if (length(keep_cells) == 0) {
  stop("No malignant cells remain after excluding rev_pathological_subtype == 'Normal'.")
}

message("Keeping ", length(keep_cells), " cells after excluding Normal pathological subtype.")
obj <- subset(obj, cells = keep_cells)
obj$final_celltype <- lineage_label
Idents(obj) <- factor(rep(lineage_label, ncol(obj)), levels = lineage_label)

message("Writing: ", output_path)
saveRDS(obj, output_path, compress = TRUE)
message("Done")
