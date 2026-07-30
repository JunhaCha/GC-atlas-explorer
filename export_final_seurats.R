#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) {
    y
  } else {
    x
  }
}

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(flag, default = NULL) {
  hit <- which(args %in% flag)
  if (length(hit) == 0 || hit[[1]] == length(args)) {
    return(default)
  }
  args[[hit[[1]] + 1L]]
}

has_flag <- function(flag) {
  any(args %in% flag)
}

file_args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(file_args) > 0) {
  sub("^--file=", "", file_args[[1]])
} else {
  file.path(getwd(), "export_final_seurats.R")
}
repo_dir <- normalizePath(dirname(script_file), winslash = "/", mustWork = FALSE)
old_wd <- getwd()
setwd(repo_dir)
on.exit(setwd(old_wd), add = TRUE)

data_dir <- normalizePath(
  arg_value("--data-dir", Sys.getenv("GC_APP_DATA_DIR", unset = dirname(repo_dir))),
  winslash = "/",
  mustWork = FALSE
)
xenium_dir <- normalizePath(
  arg_value("--xenium-dir", file.path(data_dir, "atlas_app_inputs")),
  winslash = "/",
  mustWork = FALSE
)
output_dir <- normalizePath(
  arg_value("--output-dir", file.path(data_dir, "final_seurats")),
  winslash = "/",
  mustWork = FALSE
)
overwrite <- has_flag("--overwrite")
skip_scrna <- has_flag("--skip-scrna")
skip_xenium <- has_flag("--skip-xenium")
only_scrna <- trimws(strsplit(arg_value("--scrna-lineage", ""), ",", fixed = TRUE)[[1]])
only_scrna <- only_scrna[nzchar(only_scrna)]
only_xenium <- trimws(strsplit(arg_value("--xenium-sample", ""), ",", fixed = TRUE)[[1]])
only_xenium <- only_xenium[nzchar(only_xenium)]

Sys.setenv(GC_APP_DATA_DIR = data_dir)
source(file.path(repo_dir, "R", "helpers.R"), local = TRUE)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Data directory: ", data_dir)
message("Xenium directory: ", xenium_dir)
message("Output directory: ", output_dir)

first_existing_col <- function(meta, candidates) {
  hits <- candidates[candidates %in% colnames(meta)]
  if (length(hits) == 0) {
    return(NULL)
  }
  hits[[1]]
}

is_blankish <- function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  is.na(x) | !nzchar(trimws(as.character(x))) | as.character(x) == "NA"
}

ensure_final_celltype <- function(obj, meta) {
  if ("final_celltype" %in% colnames(meta) && !all(is_blankish(meta$final_celltype))) {
    return(meta)
  }

  meta$final_celltype <- as.character(Seurat::Idents(obj))
  meta
}

derive_final_group <- function(meta) {
  patient_col <- first_existing_col(meta, c("patientID", "patient_id", "patient"))
  condition_col <- first_existing_col(meta, c("rev_condition", "condition"))
  subtype_col <- first_existing_col(
    meta,
    c("rev_pathological_subtype", "pathological_subtype", "pathological subtype")
  )

  if (is.null(patient_col) || is.null(condition_col) || is.null(subtype_col)) {
    return(rep(NA_character_, nrow(meta)))
  }

  patient <- as.character(meta[[patient_col]])
  condition <- as.character(meta[[condition_col]])
  subtype <- as.character(meta[[subtype_col]])

  primary_keep <- condition == "Primary" & !is_blankish(subtype) & subtype != "Normal"
  primary_subtype <- stats::setNames(subtype[primary_keep], patient[primary_keep])
  primary_subtype <- primary_subtype[!duplicated(names(primary_subtype))]

  revised_subtype <- subtype
  normal_keep <- condition == "Normal"
  revised_subtype[normal_keep] <- unname(primary_subtype[patient[normal_keep]])
  revised_subtype[is_blankish(revised_subtype)] <- NA_character_

  group <- paste(revised_subtype, condition, sep = "_")
  group[is.na(revised_subtype) | is_blankish(condition)] <- NA_character_
  group
}

ensure_final_group <- function(meta) {
  if ("final_group" %in% colnames(meta) && !all(is_blankish(meta$final_group))) {
    return(meta)
  }

  meta$final_group <- derive_final_group(meta)
  meta
}

curate_atlas_metadata_complete <- function(meta) {
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  out <- data.frame(row.names = rownames(meta), check.names = FALSE)

  for (spec in atlas_metadata_specs()) {
    hit <- first_existing_col(meta, spec$aliases)
    out[[spec$label]] <- if (is.null(hit)) {
      rep(NA, nrow(meta))
    } else {
      meta[[hit]]
    }
  }

  out
}

safe_output_path <- function(prefix, output_name) {
  file.path(output_dir, paste0(prefix, "_", output_name, "_full_organized.rds"))
}

write_manifest <- function(records) {
  if (length(records) == 0) {
    return(invisible(NULL))
  }

  manifest <- do.call(rbind, records)
  utils::write.table(
    manifest,
    file = file.path(output_dir, "final_seurats_manifest.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

export_scRNA_objects <- function() {
  scRNA_inputs <- list(
    list(lineage = "Entire Atlas", output_name = "Entire_Atlas", file = "seurat_merged_TME_malignant_final_umap.rds"),
    list(lineage = "Normal Epithelial", output_name = "Normal_Epithelial", file = "seurat_epithelial_normal_final_final.rds"),
    list(lineage = "Malignant", output_name = "Malignant", file = "seurat_cancercells_final.rds"),
    list(lineage = "Stromal", output_name = "Stromal", file = "seurat_Stromal_final.rds"),
    list(lineage = "CD8 T Cells", output_name = "CD8_T_Cells", file = "seurat_CD8T_final2.rds"),
    list(lineage = "CD4 T Cells", output_name = "CD4_T_Cells", file = "seurat_CD4T_final2.rds"),
    list(lineage = "B Cells", output_name = "B_Cells", file = "seurat_B_final2.rds"),
    list(lineage = "Myeloid Cells", output_name = "Myeloid_Cells", file = "seurat_Mye_final2.rds")
  )

  records <- list()

  for (entry in scRNA_inputs) {
    if (length(only_scrna) > 0 && !(entry$lineage %in% only_scrna || entry$output_name %in% only_scrna)) {
      next
    }

    input_path <- file.path(data_dir, entry$file)
    output_path <- safe_output_path("scRNA", entry$output_name)

    if (!file.exists(input_path)) {
      warning("Skipping missing full scRNA object: ", input_path, call. = FALSE)
      next
    }
    if (file.exists(output_path) && !overwrite) {
      message("Skipping existing export: ", output_path)
      next
    }

    message("\n[scRNA] Reading ", entry$lineage, ": ", input_path)
    obj <- readRDS(input_path)
    meta <- obj@meta.data
    meta <- ensure_final_celltype(obj, meta)
    meta <- ensure_final_group(meta)
    obj@meta.data <- curate_atlas_metadata_complete(meta)

    if ("final_celltype" %in% colnames(obj@meta.data) && !all(is_blankish(obj@meta.data$final_celltype))) {
      Seurat::Idents(obj) <- obj@meta.data$final_celltype
    }

    message("[scRNA] Writing ", output_path)
    saveRDS(obj, output_path)

    records[[length(records) + 1L]] <- data.frame(
      object_type = "scRNA",
      label = entry$lineage,
      raw_sample_id = NA_character_,
      source_file = input_path,
      output_file = output_path,
      cells = ncol(obj),
      features = nrow(obj),
      metadata_columns = paste(colnames(obj@meta.data), collapse = ";"),
      reductions = paste(names(obj@reductions), collapse = ";"),
      stringsAsFactors = FALSE
    )

    rm(obj, meta)
    invisible(gc())
  }

  records
}

read_xenium_manifest_for_export <- function() {
  manifest_path <- file.path(xenium_dir, "export_manifest.tsv")
  if (!file.exists(manifest_path)) {
    stop("Xenium export manifest was not found: ", manifest_path, call. = FALSE)
  }

  manifest <- utils::read.delim(manifest_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("sample_id", "cells_file", "image_meta_rds")
  missing <- setdiff(required, colnames(manifest))
  if (length(missing) > 0) {
    stop("Xenium manifest is missing column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  manifest
}

spatial_condition_from_label <- function(label) {
  if (grepl("^Normal_", label)) {
    return("Normal")
  }
  "Primary"
}

spatial_pathology_from_label <- function(label) {
  sub("_.*$", "", label)
}

export_xenium_objects <- function() {
  manifest <- read_xenium_manifest_for_export()
  label_map <- xenium_sample_labels()
  if (length(only_xenium) > 0) {
    display_labels <- unname(label_map[manifest$sample_id])
    display_labels[is.na(display_labels) | !nzchar(display_labels)] <- manifest$sample_id[is.na(display_labels) | !nzchar(display_labels)]
    keep <- manifest$sample_id %in% only_xenium | display_labels %in% only_xenium
    manifest <- manifest[keep, , drop = FALSE]
  }

  records <- list()

  for (i in seq_len(nrow(manifest))) {
    sample_id <- manifest$sample_id[[i]]
    display_label <- unname(label_map[sample_id])
    if (length(display_label) == 0 || is.na(display_label) || !nzchar(display_label)) {
      display_label <- sample_id
    }

    output_path <- file.path(output_dir, paste0("Xenium_", display_label, "_", sample_id, "_full_organized.rds"))
    if (file.exists(output_path) && !overwrite) {
      message("Skipping existing export: ", output_path)
      next
    }

    sample_dir <- file.path(xenium_dir, sample_id)
    counts_path <- file.path(sample_dir, paste0(sample_id, "_counts_sce.rds"))
    cells_path <- file.path(sample_dir, manifest$cells_file[[i]])
    image_meta_path <- file.path(sample_dir, manifest$image_meta_rds[[i]])

    if (!file.exists(counts_path)) {
      warning("Skipping missing Xenium count matrix: ", counts_path, call. = FALSE)
      next
    }
    if (!file.exists(cells_path)) {
      warning("Skipping missing Xenium cell metadata: ", cells_path, call. = FALSE)
      next
    }
    if (!file.exists(image_meta_path)) {
      warning("Skipping missing Xenium image metadata: ", image_meta_path, call. = FALSE)
      next
    }

    message("\n[Xenium] Reading ", display_label, " (", sample_id, ")")
    counts <- readRDS(counts_path)
    cells <- read_xenium_table_file(cells_path)
    image_meta <- readRDS(image_meta_path)
    cells <- xenium_transform_cells(cells, image_meta)

    cells <- as.data.frame(cells, stringsAsFactors = FALSE)
    rownames(cells) <- cells$cell_id

    shared_cells <- intersect(colnames(counts), rownames(cells))
    if (length(shared_cells) == 0) {
      stop("No shared cell IDs between counts and cell metadata for sample ", sample_id, ".", call. = FALSE)
    }
    if (length(shared_cells) != ncol(counts) || length(shared_cells) != nrow(cells)) {
      message("[Xenium] Aligning ", length(shared_cells), " shared cells for ", sample_id)
    }

    counts <- counts[, shared_cells, drop = FALSE]
    cells <- cells[shared_cells, , drop = FALSE]

    spatial_meta <- data.frame(
      sample_id = display_label,
      raw_sample_id = sample_id,
      x = as.numeric(cells$x),
      y = as.numeric(cells$y),
      x_plot = as.numeric(cells$x_plot),
      y_plot = as.numeric(cells$y_plot),
      celltype = as.character(cells$celltype),
      subtype = as.character(cells$subtype),
      tumor_normal = as.character(cells$tumor_normal),
      direct_neighborhood = as.character(cells$neighborhood),
      condition = spatial_condition_from_label(display_label),
      `pathological subtype` = spatial_pathology_from_label(display_label),
      row.names = shared_cells,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    obj <- Seurat::CreateSeuratObject(
      counts = counts,
      assay = "Xenium",
      project = display_label,
      meta.data = spatial_meta
    )

    Seurat::Idents(obj) <- obj@meta.data$celltype
    embeddings <- as.matrix(obj@meta.data[, c("x", "y"), drop = FALSE])
    colnames(embeddings) <- c("spatial_1", "spatial_2")
    obj[["spatial"]] <- Seurat::CreateDimReducObject(
      embeddings = embeddings,
      key = "spatial_",
      assay = Seurat::DefaultAssay(obj)
    )

    message("[Xenium] Writing ", output_path)
    saveRDS(obj, output_path)

    records[[length(records) + 1L]] <- data.frame(
      object_type = "Xenium",
      label = display_label,
      raw_sample_id = sample_id,
      source_file = counts_path,
      output_file = output_path,
      cells = ncol(obj),
      features = nrow(obj),
      metadata_columns = paste(colnames(obj@meta.data), collapse = ";"),
      reductions = paste(names(obj@reductions), collapse = ";"),
      stringsAsFactors = FALSE
    )

    rm(obj, counts, cells, spatial_meta, embeddings)
    invisible(gc())
  }

  records
}

all_records <- list()
if (!skip_scrna) {
  all_records <- c(all_records, export_scRNA_objects())
}
if (!skip_xenium) {
  all_records <- c(all_records, export_xenium_objects())
}

write_manifest(all_records)
message("\nFinished final Seurat exports.")
