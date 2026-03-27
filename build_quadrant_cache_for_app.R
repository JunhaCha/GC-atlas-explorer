#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg_value <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) {
    return(default)
  }
  args[[idx + 1]]
}

input_path <- get_arg_value("--input")
output_path <- get_arg_value("--output")
slot_name <- get_arg_value("--slot", "data")
celltype_col <- get_arg_value("--celltype-col", "final_celltype")
group_col <- get_arg_value("--group-col", "final_group")

if (is.null(input_path) || is.null(output_path)) {
  stop("Usage: Rscript build_quadrant_cache_for_app.R --input input.rds --output output_quadrant_cache.rds [--slot data] [--celltype-col final_celltype] [--group-col final_group]")
}

source("/Users/junhacha/Documents/Playground/seurat-shiny-explorer/R/helpers.R")

message("Reading Seurat object: ", input_path)
obj <- readRDS(input_path)
if (!inherits(obj, "Seurat")) {
  stop("Input file is not a Seurat object: ", input_path)
}

message("Building quadrant FC cache")
cache <- build_quadrant_fc_cache(
  obj = obj,
  genes_to_use = rownames(obj),
  assay_name = DefaultAssay(obj),
  slot_name = slot_name,
  celltype_col = celltype_col,
  group_col = group_col
)

message("Saving cache: ", output_path)
saveRDS(cache, output_path, compress = TRUE)
message("Done")
