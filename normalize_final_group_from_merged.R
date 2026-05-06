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

reference_path <- if (length(args) >= 2 && nzchar(args[[2]])) {
  args[[2]]
} else {
  file.path(data_dir, "seurat_merged_TME_malignant_final_umap_app_slim.rds")
}

target_files <- if (length(args) >= 3) {
  normalizePath(args[3:length(args)], winslash = "/", mustWork = FALSE)
} else {
  character(0)
}

group_levels <- c(
  "Diffuse_Primary",
  "Diffuse_Normal",
  "Intestinal_Primary",
  "Intestinal_Normal",
  "Mixed_Primary",
  "Mixed_Normal"
)
subtype_levels <- c("Diffuse", "Intestinal", "Mixed")

if (!file.exists(reference_path)) {
  stop("Reference merged Seurat object not found: ", reference_path)
}

message("Reading reference object: ", reference_path)
reference_obj <- readRDS(reference_path)
if (!inherits(reference_obj, "Seurat")) {
  stop("Reference file is not a Seurat object: ", reference_path)
}

ref_md <- reference_obj@meta.data
required_cols <- c("patientID", "rev_condition", "rev_pathological_subtype")
missing_cols <- setdiff(required_cols, colnames(ref_md))
if (length(missing_cols) > 0) {
  stop("Reference object is missing required columns: ", paste(missing_cols, collapse = ", "))
}

primary_map <- unique(ref_md[ref_md$rev_condition == "Primary", c("patientID", "rev_pathological_subtype"), drop = FALSE])
primary_map <- primary_map[
  !is.na(primary_map$patientID) &
    nzchar(as.character(primary_map$patientID)) &
    !is.na(primary_map$rev_pathological_subtype) &
    primary_map$rev_pathological_subtype %in% subtype_levels,
  ,
  drop = FALSE
]

patient_to_subtype <- split(as.character(primary_map$rev_pathological_subtype), as.character(primary_map$patientID))
patient_to_subtype <- lapply(patient_to_subtype, unique)
conflicts <- names(patient_to_subtype)[vapply(patient_to_subtype, length, integer(1)) > 1]
if (length(conflicts) > 0) {
  stop(
    "Some patients map to multiple primary subtypes in the reference object. Examples: ",
    paste(head(conflicts, 10), collapse = ", ")
  )
}
patient_to_subtype <- vapply(patient_to_subtype, `[`, character(1), 1)

files <- if (length(target_files) > 0) {
  target_files
} else {
  files_all <- sort(Sys.glob(file.path(data_dir, "seurat*.rds")))
  files_all[!(endsWith(files_all, "_avg_cache.rds") | endsWith(files_all, "_quadrant_cache.rds"))]
}

updated <- character(0)

for (f in files) {
  message("")
  message("[normalize-final-group] Checking ", basename(f))
  obj <- readRDS(f)

  if (!inherits(obj, "Seurat")) {
    message("  Not a Seurat object, skipping.")
    next
  }

  md <- obj@meta.data
  missing_cols <- setdiff(required_cols, colnames(md))
  if (length(missing_cols) > 0) {
    message("  Missing required columns: ", paste(missing_cols, collapse = ", "), ". Skipping.")
    next
  }

  primary_subtype <- unname(patient_to_subtype[as.character(md$patientID)])
  revised_subtype <- ifelse(md$rev_condition == "Normal", primary_subtype, as.character(md$rev_pathological_subtype))
  revised_subtype[!(revised_subtype %in% subtype_levels)] <- NA_character_

  keep_cells <- !is.na(revised_subtype) & md$rev_condition %in% c("Primary", "Normal")
  dropped_cells <- sum(!keep_cells)
  if (dropped_cells > 0) {
    message("  Dropping ", dropped_cells, " cells without mapped Diffuse/Intestinal/Mixed subtype.")
  }

  if (!all(keep_cells)) {
    obj <- subset(obj, cells = rownames(md)[keep_cells])
    md <- obj@meta.data
    primary_subtype <- primary_subtype[keep_cells]
    revised_subtype <- revised_subtype[keep_cells]
  }

  final_group <- paste(revised_subtype, as.character(md$rev_condition), sep = "_")
  md$final_group <- factor(final_group, levels = group_levels)
  obj@meta.data <- md

  saveRDS(obj, f, compress = TRUE)
  message("  Saved normalized final_group to ", basename(f))
  updated <- c(updated, basename(f))
}

message("")
message("UPDATED")
if (length(updated) > 0) {
  message(paste(updated, collapse = "\n"))
} else {
  message("none")
}
