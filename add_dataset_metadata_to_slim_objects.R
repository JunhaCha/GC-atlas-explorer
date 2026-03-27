#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

slim_paths <- c(
  "/Users/junhacha/Documents/Playground/seurat_merged_TME_malignant_final_umap_app_slim.rds",
  "/Users/junhacha/Documents/Playground/seurat_epithelial_normal_final_final_app_slim.rds",
  "/Users/junhacha/Documents/Playground/seurat_cancercells_final_app_slim.rds",
  "/Users/junhacha/Documents/Playground/seurat_Stromal_final_app_slim.rds",
  "/Users/junhacha/Documents/Playground/seurat_CD8T_final2_app_slim.rds",
  "/Users/junhacha/Documents/Playground/seurat_CD4T_final2_app_slim.rds",
  "/Users/junhacha/Documents/Playground/seurat_B_final2_app_slim.rds",
  "/Users/junhacha/Documents/Playground/seurat_Mye_final2_app_slim.rds"
)

dataset_levels_requested <- c(
  "Sathe_2020",
  "Jeong_2021",
  "Kang_2022",
  "Sun_2022",
  "Wang_2023",
  "Sathe_2023",
  "Cheng_2024"
)

parse_dataset_from_orig_ident <- function(orig_ident) {
  out <- sub(":.*$", "", as.character(orig_ident))
  out[is.na(orig_ident) | !nzchar(as.character(orig_ident))] <- NA_character_
  out
}

infer_dataset_from_sample_id <- function(sample_id) {
  x <- as.character(sample_id)
  out <- rep(NA_character_, length(x))

  out[grepl("^GSM510", x)] <- "Jeong_2021"
  out[grepl("^[0-9]{6}[NT]([12])?$", x)] <- "Kang_2022"
  out[grepl("^GC[0-9]+-", x)] <- "Sun_2022"
  out[grepl("^GSM45463", x)] <- "Kim_2022"
  out[grepl("^P[0-9]+-(Ad|Ca)$", x)] <- "Wang_2023"
  out[grepl("^Pt[0-9]+_(Ad|P)$", x)] <- "Cheng_2024"
  out[grepl("^sample[0-9]+$", x)] <- "Kumar_2022"
  out[grepl("^P[0-9]+-(N|T|N1|N2|T1|T2)$", x)] <- "Sathe_2020"
  out[grepl("^P8654-T$", x)] <- "Sathe_2023"

  out
}

sample_lookup <- c()

for (path in slim_paths) {
  if (!file.exists(path)) {
    next
  }

  obj <- readRDS(path)
  md <- obj@meta.data

  if ("orig.ident" %in% colnames(md) && "sampleID" %in% colnames(md)) {
    dataset_vals <- parse_dataset_from_orig_ident(md$orig.ident)
    keep <- !is.na(dataset_vals) & nzchar(dataset_vals) &
      !is.na(md$sampleID) & nzchar(as.character(md$sampleID))
    if (any(keep)) {
      sample_lookup <- c(
        sample_lookup,
        stats::setNames(as.character(dataset_vals[keep]), as.character(md$sampleID[keep]))
      )
    }
  }

  if ("dataset" %in% colnames(md) && "sampleID" %in% colnames(md)) {
    keep <- !is.na(md$dataset) & nzchar(as.character(md$dataset)) &
      !is.na(md$sampleID) & nzchar(as.character(md$sampleID))
    if (any(keep)) {
      sample_lookup <- c(
        sample_lookup,
        stats::setNames(as.character(md$dataset[keep]), as.character(md$sampleID[keep]))
      )
    }
  }
}

sample_lookup <- sample_lookup[!duplicated(names(sample_lookup))]

for (path in slim_paths) {
  if (!file.exists(path)) {
    message("Skipping missing slim object: ", path)
    next
  }

  message("Updating dataset metadata in: ", basename(path))
  obj <- readRDS(path)
  md <- obj@meta.data

  dataset_vals <- NULL
  if ("orig.ident" %in% colnames(md)) {
    dataset_vals <- parse_dataset_from_orig_ident(md$orig.ident)
  } else if ("sampleID" %in% colnames(md)) {
    dataset_vals <- unname(sample_lookup[as.character(md$sampleID)])
    inferred_vals <- infer_dataset_from_sample_id(md$sampleID)
    missing_lookup <- is.na(dataset_vals) | !nzchar(dataset_vals)
    dataset_vals[missing_lookup] <- inferred_vals[missing_lookup]
  } else {
    stop("Cannot derive dataset for ", basename(path), ": missing both orig.ident and sampleID.")
  }

  missing_dataset <- is.na(dataset_vals) | !nzchar(dataset_vals)
  if (any(missing_dataset)) {
    unresolved_samples <- unique(as.character(md$sampleID[missing_dataset]))
    unresolved_samples <- unresolved_samples[!is.na(unresolved_samples) & nzchar(unresolved_samples)]
    stop(
      "Could not resolve dataset for some cells in ", basename(path),
      if (length(unresolved_samples) > 0) paste0(". Example sampleID values: ", paste(head(unresolved_samples, 10), collapse = ", ")) else "."
    )
  }

  extra_levels <- sort(setdiff(unique(dataset_vals), dataset_levels_requested))
  md$dataset <- factor(dataset_vals, levels = c(dataset_levels_requested, extra_levels))
  obj@meta.data <- md

  saveRDS(obj, path, compress = TRUE)
}

message("Finished updating dataset metadata in slim objects.")
