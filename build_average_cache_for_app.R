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

pick_sample_column <- function(md) {
  candidates <- c("sampleID", "sample", "orig.ident", "sample_id")
  hits <- candidates[candidates %in% colnames(md)]
  if (length(hits) == 0) {
    stop("Could not find a sample column. Expected one of: ", paste(candidates, collapse = ", "))
  }
  hits[[1]]
}

derive_pathological_group <- function(md) {
  if ("rev_pathological_group" %in% colnames(md)) {
    return(as.character(md$rev_pathological_group))
  }
  if ("final_group" %in% colnames(md)) {
    x <- as.character(md$final_group)
    out <- rep(NA_character_, length(x))
    out[grepl("_Normal$", x)] <- "Normal"
    out[grepl("^Diffuse", x)] <- "Diffuse"
    out[grepl("^Intestinal", x)] <- "Intestinal"
    out[grepl("^Mixed", x)] <- "Mixed"
    return(out)
  }
  rep(NA_character_, nrow(md))
}

mode_non_na <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) {
    return(NA_character_)
  }
  tab <- sort(table(x), decreasing = TRUE)
  names(tab)[[1]]
}

args <- parse_args()
input_path <- get_arg(args, "input", required = TRUE)
output_path <- get_arg(args, "output", required = TRUE)
assay_name <- get_arg(args, "assay", default = NULL)
layer_name <- get_arg(args, "layer", default = "data")

message("Reading: ", input_path)
obj <- readRDS(input_path)
if (!inherits(obj, "Seurat")) {
  stop("Input file is not a Seurat object.")
}

assay_to_use <- assay_name %||% DefaultAssay(obj)
sample_col <- pick_sample_column(obj@meta.data)
celltype_col <- if ("final_celltype" %in% colnames(obj@meta.data)) "final_celltype" else NULL

md <- obj@meta.data
md$rev_pathological_group <- derive_pathological_group(md)
if (is.null(celltype_col)) {
  md$final_celltype <- as.character(Idents(obj))
  celltype_col <- "final_celltype"
}

keep_cells <- !is.na(md[[sample_col]]) & nzchar(as.character(md[[sample_col]])) &
  !is.na(md[[celltype_col]]) & nzchar(as.character(md[[celltype_col]]))
md <- md[keep_cells, , drop = FALSE]
obj <- subset(obj, cells = rownames(md))
DefaultAssay(obj) <- assay_to_use

group_key <- paste(md[[sample_col]], md[[celltype_col]], sep = " || ")
group_levels <- unique(group_key)
pb_ids <- sprintf("pb%06d", seq_along(group_levels))
pb_map <- stats::setNames(pb_ids, group_levels)
obj$pseudobulk_id <- unname(pb_map[group_key])
avg_mat <- AverageExpression(
  object = obj,
  assays = assay_to_use,
  group.by = "pseudobulk_id",
  layer = layer_name,
  verbose = FALSE
)[[assay_to_use]]
avg_mat <- as.matrix(avg_mat)

metadata_fields <- c(
  sample_col,
  celltype_col,
  "dataset",
  "rev_pathological_group",
  "rev_pathological_subtype",
  "rev_condition",
  "rev_molecular_subtype",
  "patientID",
  "TNM",
  "Phase",
  "final_group"
)
metadata_fields <- unique(metadata_fields[metadata_fields %in% colnames(md)])

meta_out <- data.frame(
  pseudobulk_id = colnames(avg_mat),
  stringsAsFactors = FALSE
)
inverse_pb_map <- stats::setNames(names(pb_map), unname(pb_map))
split_keys <- strsplit(unname(inverse_pb_map[colnames(avg_mat)]), " \\|\\| ", perl = TRUE)
meta_out[[sample_col]] <- vapply(split_keys, `[`, "", 1)
meta_out[[celltype_col]] <- vapply(split_keys, `[`, "", 2)

for (field in setdiff(metadata_fields, c(sample_col, celltype_col))) {
  vals <- vapply(colnames(avg_mat), function(pb_id) {
    source_key <- inverse_pb_map[[pb_id]]
    mode_non_na(md[group_key == source_key, field, drop = TRUE])
  }, character(1))
  meta_out[[field]] <- vals
}

if (!"sampleID" %in% colnames(meta_out) && sample_col != "sampleID") {
  meta_out$sampleID <- meta_out[[sample_col]]
}
if (!"final_celltype" %in% colnames(meta_out) && celltype_col != "final_celltype") {
  meta_out$final_celltype <- meta_out[[celltype_col]]
}

rownames(meta_out) <- meta_out$pseudobulk_id
ord_fields <- c("final_celltype", "rev_pathological_group", "dataset", "sampleID", "patientID", "final_group", "pseudobulk_id")
ord_fields <- ord_fields[ord_fields %in% colnames(meta_out)]
ord <- do.call(order, unname(lapply(ord_fields, function(field) as.character(meta_out[[field]]))))
meta_out <- meta_out[ord, , drop = FALSE]
avg_mat <- avg_mat[, rownames(meta_out), drop = FALSE]

cache <- list(
  avg_mat = avg_mat,
  meta = meta_out,
  assay = assay_to_use,
  layer = layer_name,
  source_path = input_path,
  sample_col = sample_col,
  celltype_col = celltype_col
)

message("Writing: ", output_path)
saveRDS(cache, output_path, compress = TRUE)
message(sprintf("Cache dimensions: %d genes x %d pseudobulk columns", nrow(avg_mat), ncol(avg_mat)))
