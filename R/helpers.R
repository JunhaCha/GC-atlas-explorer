suppressPackageStartupMessages({
  library(shiny)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

app_code_dir <- function() {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

app_data_dir <- function() {
  env_dir <- Sys.getenv("GC_APP_DATA_DIR", unset = "")
  if (nzchar(env_dir)) {
    return(normalizePath(env_dir, winslash = "/", mustWork = FALSE))
  }

  local_data_dir <- file.path(app_code_dir(), "data")
  if (dir.exists(local_data_dir)) {
    return(normalizePath(local_data_dir, winslash = "/", mustWork = FALSE))
  }

  normalizePath(dirname(app_code_dir()), winslash = "/", mustWork = FALSE)
}

app_data_path <- function(...) {
  normalizePath(file.path(app_data_dir(), ...), winslash = "/", mustWork = FALSE)
}

app_script_path <- function(script_name) {
  normalizePath(file.path(app_code_dir(), script_name), winslash = "/", mustWork = FALSE)
}

atlas_object_choices_default <- function() {
  c(
    "Entire Atlas" = prefer_app_slim_path(app_data_path("seurat_merged_TME_malignant_final_umap.rds")),
    "Normal Epithelial" = prefer_app_slim_path(app_data_path("seurat_epithelial_normal_final_final.rds")),
    "Malignant" = prefer_app_slim_path(app_data_path("seurat_cancercells_final.rds")),
    "Stromal" = prefer_app_slim_path(app_data_path("seurat_Stromal_final.rds")),
    "CD8 T Cells" = prefer_app_slim_path(app_data_path("seurat_CD8T_final2.rds")),
    "CD4 T Cells" = prefer_app_slim_path(app_data_path("seurat_CD4T_final2.rds")),
    "B Cells" = prefer_app_slim_path(app_data_path("seurat_B_final2.rds")),
    "Myeloid Cells" = prefer_app_slim_path(app_data_path("seurat_Mye_final2.rds"))
  )
}

xenium_object_choices_default <- function() {
  c("Xenium object not configured yet" = "")
}

xenium_app_inputs_dir <- function() {
  candidates <- c(
    app_data_path("atlas_app_inputs"),
    app_data_dir()
  )
  candidates <- unique(candidates[file.exists(candidates) | dir.exists(candidates)])

  manifest_hits <- candidates[file.exists(file.path(candidates, "export_manifest.tsv"))]
  if (length(manifest_hits) > 0) {
    return(normalizePath(manifest_hits[[1]], winslash = "/", mustWork = FALSE))
  }

  normalizePath(candidates[[1]] %||% app_data_path("atlas_app_inputs"), winslash = "/", mustWork = FALSE)
}

xenium_manifest_path <- function() {
  normalizePath(file.path(xenium_app_inputs_dir(), "export_manifest.tsv"), winslash = "/", mustWork = FALSE)
}

xenium_curated_gene_list_path <- function() {
  normalizePath(
    file.path(app_code_dir(), "config", "xenium_curated_gene_list.csv"),
    winslash = "/",
    mustWork = FALSE
  )
}

xenium_sample_labels_path <- function() {
  normalizePath(
    file.path(app_code_dir(), "config", "xenium_sample_labels.csv"),
    winslash = "/",
    mustWork = FALSE
  )
}

xenium_sample_labels <- function() {
  path <- xenium_sample_labels_path()
  if (!file.exists(path)) {
    return(stats::setNames(character(0), character(0)))
  }

  labels_df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c("sample_id", "display_label")
  if (!all(required_cols %in% colnames(labels_df))) {
    stop(
      paste0(
        "Xenium sample-label file must contain columns: ",
        paste(required_cols, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  keep <- !is.na(labels_df$sample_id) & nzchar(labels_df$sample_id) &
    !is.na(labels_df$display_label) & nzchar(labels_df$display_label)
  labels_df <- labels_df[keep, required_cols, drop = FALSE]

  stats::setNames(labels_df$display_label, labels_df$sample_id)
}

xenium_curated_gene_table <- function() {
  path <- xenium_curated_gene_list_path()
  if (!file.exists(path)) {
    stop(paste0("Xenium curated gene list not found at ", path, "."), call. = FALSE)
  }

  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

xenium_curated_genes <- function() {
  unique(xenium_curated_gene_table()$gene)
}

xenium_counts_matrix_path <- function(sample_id) {
  normalizePath(
    file.path(xenium_app_inputs_dir(), sample_id, paste0(sample_id, "_counts_sce.rds")),
    winslash = "/",
    mustWork = FALSE
  )
}

xenium_sample_manifest <- function() {
  manifest_path <- xenium_manifest_path()
  if (!file.exists(manifest_path)) {
    stop(paste0("Xenium manifest not found at ", manifest_path, "."), call. = FALSE)
  }

  manifest <- utils::read.delim(
    manifest_path,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!("sample_id" %in% colnames(manifest))) {
    stop("Xenium manifest must contain a 'sample_id' column.", call. = FALSE)
  }

  path_fields <- c(
    "cells_file",
    "image_file",
    "image_meta_json",
    "image_meta_rds",
    "neighborhood_summary_file",
    "neighborhood_composition_file"
  )

  root_dir <- dirname(manifest_path)
  for (field in intersect(path_fields, colnames(manifest))) {
    manifest[[field]] <- normalizePath(
      file.path(root_dir, manifest$sample_id, manifest[[field]]),
      winslash = "/",
      mustWork = FALSE
    )
  }

  manifest
}

xenium_sample_choices <- function(manifest = xenium_sample_manifest()) {
  sample_ids <- as.character(manifest$sample_id)
  label_map <- xenium_sample_labels()
  display_labels <- unname(label_map[sample_ids])
  display_labels[is.na(display_labels) | !nzchar(display_labels)] <- sample_ids[is.na(display_labels) | !nzchar(display_labels)]

  parse_group_rank <- function(labels, pattern, rank) {
    hits <- grepl(pattern, labels)
    idx <- rep(Inf, length(labels))
    if (any(hits)) {
      idx[hits] <- suppressWarnings(as.numeric(sub(pattern, "\\1", labels[hits])))
    }
    list(group = ifelse(hits, rank, Inf), order = idx)
  }

  diffuse <- parse_group_rank(display_labels, "^Diffuse_(\\d+)$", 1)
  intestinal <- parse_group_rank(display_labels, "^Intestinal_(\\d+)$", 2)
  normal <- parse_group_rank(display_labels, "^Normal_(\\d+)$", 3)

  group_rank <- pmin(diffuse$group, intestinal$group, normal$group)
  within_rank <- pmin(diffuse$order, intestinal$order, normal$order)
  fallback_label <- tolower(display_labels)
  ord <- order(group_rank, within_rank, fallback_label, sample_ids, na.last = TRUE)

  stats::setNames(sample_ids[ord], display_labels[ord])
}

xenium_manifest_row <- function(sample_id, manifest = xenium_sample_manifest()) {
  if (is.null(sample_id) || !nzchar(sample_id)) {
    stop("Choose a Xenium sample.", call. = FALSE)
  }
  row <- manifest[manifest$sample_id %in% sample_id, , drop = FALSE]
  if (nrow(row) != 1) {
    stop(paste0("Sample '", sample_id, "' was not found in the Xenium manifest."), call. = FALSE)
  }
  row[1, , drop = FALSE]
}

xenium_image_meta <- function(sample_row) {
  meta_path <- sample_row$image_meta_rds[[1]] %||% ""
  if (!file.exists(meta_path)) {
    stop(paste0("Image metadata file was not found for sample ", sample_row$sample_id[[1]], "."), call. = FALSE)
  }
  readRDS(meta_path)
}

read_xenium_table_file <- function(path) {
  ext <- tolower(tools::file_ext(path %||% ""))

  fallback_path <- function(new_ext) {
    sub(paste0("\\.", ext, "$"), paste0(".", new_ext), path, ignore.case = TRUE)
  }

  if (identical(ext, "rds")) {
    return(readRDS(path))
  }

  if (identical(ext, "qs")) {
    if (requireNamespace("qs", quietly = TRUE)) {
      return(qs::qread(path))
    }

    alt_rds <- fallback_path("rds")
    if (file.exists(alt_rds)) {
      return(readRDS(alt_rds))
    }

    stop("Install the 'qs' package to read Xenium .qs tables, or provide matching .rds fallbacks.", call. = FALSE)
  }

  if (identical(ext, "parquet")) {
    if (requireNamespace("arrow", quietly = TRUE)) {
      return(as.data.frame(arrow::read_parquet(path), stringsAsFactors = FALSE))
    }

    alt_rds <- fallback_path("rds")
    if (file.exists(alt_rds)) {
      return(readRDS(alt_rds))
    }

    stop("Install the 'arrow' package to read Xenium .parquet tables, or provide matching .rds fallbacks.", call. = FALSE)
  }

  stop(paste0("Unsupported Xenium table format: ", path), call. = FALSE)
}

xenium_transform_cells <- function(cells_df, image_meta) {
  if (!all(c("x", "y") %in% colnames(cells_df))) {
    stop("Xenium cell table must contain x and y columns.", call. = FALSE)
  }

  x_range <- as.numeric(image_meta$x_range %||% c(min(cells_df$x, na.rm = TRUE), max(cells_df$x, na.rm = TRUE)))
  y_range <- as.numeric(image_meta$y_range %||% c(min(cells_df$y, na.rm = TRUE), max(cells_df$y, na.rm = TRUE)))
  display_width <- as.numeric(image_meta$displayed_width %||% image_meta$original_width %||% 1)
  display_height <- as.numeric(image_meta$displayed_height %||% image_meta$original_height %||% 1)

  x_span <- diff(x_range)
  y_span <- diff(y_range)
  if (!is.finite(x_span) || x_span == 0) {
    x_span <- 1
  }
  if (!is.finite(y_span) || y_span == 0) {
    y_span <- 1
  }

  x_scaled <- (as.numeric(cells_df$x) - x_range[[1]]) / x_span * display_width
  y_scaled <- (as.numeric(cells_df$y) - y_range[[1]]) / y_span * display_height

  if (isTRUE(image_meta$flip_y %||% FALSE)) {
    y_scaled <- display_height - y_scaled
  }

  dplyr::mutate(
    cells_df,
    x_plot = pmin(pmax(x_scaled, 0), display_width),
    y_plot = pmin(pmax(y_scaled, 0), display_height)
  )
}

xenium_sample_bundle <- function(sample_id, manifest = xenium_sample_manifest()) {
  sample_row <- xenium_manifest_row(sample_id, manifest)

  cells_path <- sample_row$cells_file[[1]] %||% ""
  image_path <- sample_row$image_file[[1]] %||% ""
  summary_path <- sample_row$neighborhood_summary_file[[1]] %||% ""
  composition_path <- sample_row$neighborhood_composition_file[[1]] %||% ""

  if (!file.exists(cells_path)) {
    stop(paste0("Cell table was not found for sample ", sample_id, "."), call. = FALSE)
  }
  if (!file.exists(image_path)) {
    stop(paste0("Background image was not found for sample ", sample_id, "."), call. = FALSE)
  }
  if (!file.exists(summary_path)) {
    stop(paste0("Neighborhood summary was not found for sample ", sample_id, "."), call. = FALSE)
  }
  if (!file.exists(composition_path)) {
    stop(paste0("Neighborhood composition was not found for sample ", sample_id, "."), call. = FALSE)
  }

  image_meta <- xenium_image_meta(sample_row)
  cells_df <- xenium_transform_cells(read_xenium_table_file(cells_path), image_meta)

  list(
    sample_id = sample_id,
    manifest_row = sample_row,
    cells = cells_df,
    image_path = image_path,
    image_meta = image_meta,
    neighborhood_summary = read_xenium_table_file(summary_path),
    neighborhood_composition = read_xenium_table_file(composition_path)
  )
}

xenium_marker_dir <- function(sample_id) {
  normalizePath(
    file.path(xenium_app_inputs_dir(), sample_id, "marker_expr"),
    winslash = "/",
    mustWork = FALSE
  )
}

xenium_marker_manifest_path <- function(sample_id) {
  normalizePath(
    file.path(xenium_marker_dir(sample_id), "marker_genes.rds"),
    winslash = "/",
    mustWork = FALSE
  )
}

xenium_marker_manifest <- function(sample_id) {
  path <- xenium_marker_manifest_path(sample_id)
  if (!file.exists(path)) {
    stop(
      paste0(
        "Marker assets were not found for sample ",
        sample_id,
        ". Run ",
        app_script_path("build_xenium_marker_assets.R"),
        " first."
      ),
      call. = FALSE
    )
  }
  readRDS(path)
}

xenium_marker_available_genes <- function(sample_id) {
  xenium_marker_manifest(sample_id)$gene
}

xenium_marker_expression_path <- function(sample_id, gene) {
  manifest <- xenium_marker_manifest(sample_id)
  row <- manifest[manifest$gene %in% gene, , drop = FALSE]
  if (nrow(row) != 1) {
    stop(paste0("Gene '", gene, "' is not available for sample ", sample_id, "."), call. = FALSE)
  }
  normalizePath(file.path(xenium_marker_dir(sample_id), row$file_name[[1]]), winslash = "/", mustWork = FALSE)
}

xenium_read_marker_expression <- function(bundle, gene) {
  expr_path <- xenium_marker_expression_path(bundle$sample_id, gene)
  marker_obj <- readRDS(expr_path)

  n_cells <- nrow(bundle$cells)
  values <- numeric(n_cells)

  if (length(marker_obj$index) > 0) {
    values[marker_obj$index] <- marker_obj$value
  }

  values
}

color_rds_paths <- c(
  dataset = app_data_path("dataset_group_colors.rds"),
  path_group = app_data_path("path_group_colors.rds"),
  celltype = app_data_path("subcluster_ct_colors_combined.rds")
)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) {
    y
  } else {
    x
  }
}

safe_seurat_read <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  bundle_path <- matching_atlas_bundle_path(path)
  validate(
    need(
      !identical(bundle_path, path) && file.exists(bundle_path),
      paste0(
        "Atlas runtime bundle was not found for this lineage. Build it first:\n",
        "Rscript ", app_script_path("build_all_atlas_runtime_bundles.R")
      )
    )
  )
  cache_key <- bundle_path
  if (exists(cache_key, envir = .gc_seurat_object_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .gc_seurat_object_cache, inherits = FALSE))
  }
  obj <- cached_read_rds(bundle_path)
  validate(need(is_atlas_runtime_bundle(obj), "Atlas runtime bundle file is not in the expected format."))
  assign(cache_key, obj, envir = .gc_seurat_object_cache)
  obj
}

.gc_seurat_object_cache <- new.env(parent = emptyenv())
.gc_rds_cache <- new.env(parent = emptyenv())
.gc_parquet_cache <- new.env(parent = emptyenv())
.gc_feature_cache <- new.env(parent = emptyenv())

matching_atlas_bundle_path <- function(path) {
  if (is.null(path) || !nzchar(path)) {
    return(path)
  }
  v2_path <- sub("\\.rds$", "_atlas_bundle_v2.rds", path, ignore.case = TRUE)
  if (file.exists(v2_path)) {
    return(v2_path)
  }
  sub("\\.rds$", "_atlas_bundle.rds", path, ignore.case = TRUE)
}

is_atlas_runtime_bundle <- function(obj) {
  inherits(obj, "atlas_runtime_bundle")
}

atlas_bundle_manifest <- function(obj) {
  unclass(obj)
}

cached_read_parquet <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (exists(path, envir = .gc_parquet_cache, inherits = FALSE)) {
    return(get(path, envir = .gc_parquet_cache, inherits = FALSE))
  }
  validate(need(requireNamespace("arrow", quietly = TRUE), "Install the 'arrow' package to read atlas parquet runtime files."))
  value <- as.data.frame(arrow::read_parquet(path), stringsAsFactors = FALSE)
  assign(path, value, envir = .gc_parquet_cache)
  value
}

atlas_bundle_meta_df <- function(obj) {
  manifest <- atlas_bundle_manifest(obj)
  if (!is.null(manifest[["meta_override"]])) {
    return(manifest[["meta_override"]])
  }
  if (!is.null(manifest[["meta"]])) {
    return(as.data.frame(manifest[["meta"]], stringsAsFactors = FALSE))
  }

  meta_path <- manifest[["meta_path"]] %||% NULL
  validate(need(!is.null(meta_path) && file.exists(meta_path), "Metadata parquet for the atlas runtime bundle was not found."))
  meta <- cached_read_parquet(meta_path)
  validate(need("cell_id" %in% colnames(meta), "Atlas metadata parquet must contain a cell_id column."))
  rownames(meta) <- as.character(meta$cell_id)
  meta$cell_id <- NULL
  meta
}

atlas_bundle_base_cells <- function(obj) {
  manifest <- atlas_bundle_manifest(obj)
  meta_path <- manifest[["meta_path"]] %||% NULL
  validate(need(!is.null(meta_path) && file.exists(meta_path), "Metadata parquet for the atlas runtime bundle was not found."))
  meta <- cached_read_parquet(meta_path)
  validate(need("cell_id" %in% colnames(meta), "Atlas metadata parquet must contain a cell_id column."))
  as.character(meta$cell_id)
}

atlas_bundle_umap_df <- function(obj) {
  manifest <- atlas_bundle_manifest(obj)
  if (!is.null(manifest[["umap_override"]])) {
    return(manifest[["umap_override"]])
  }
  if (!is.null(manifest[["umap"]])) {
    coords <- as.data.frame(manifest[["umap"]], stringsAsFactors = FALSE)
    validate(need(ncol(coords) >= 2, "The atlas runtime bundle does not contain two-dimensional UMAP coordinates."))
    colnames(coords)[1:2] <- c("UMAP_1", "UMAP_2")
    return(coords[, 1:2, drop = FALSE])
  }

  umap_path <- manifest[["umap_path"]] %||% NULL
  validate(need(!is.null(umap_path) && file.exists(umap_path), "UMAP parquet for the atlas runtime bundle was not found."))
  coords <- cached_read_parquet(umap_path)
  validate(need(all(c("cell_id", "UMAP_1", "UMAP_2") %in% colnames(coords)), "Atlas UMAP parquet must contain cell_id, UMAP_1, and UMAP_2 columns."))
  rownames(coords) <- as.character(coords$cell_id)
  coords <- coords[, c("UMAP_1", "UMAP_2"), drop = FALSE]
  coords
}

atlas_meta <- function(obj) {
  if (is_atlas_runtime_bundle(obj)) {
    return(atlas_bundle_meta_df(obj))
  }
  obj[[]]
}

atlas_default_assay <- function(obj) {
  if (is_atlas_runtime_bundle(obj)) {
    return(atlas_bundle_manifest(obj)[["assay"]] %||% "RNA")
  }
  DefaultAssay(obj)
}

atlas_reduction_name <- function(obj) {
  if (is_atlas_runtime_bundle(obj)) {
    return(atlas_bundle_manifest(obj)[["reduction_name"]] %||% "UMAP")
  }
  preferred_umap_reduction(obj)
}

atlas_features <- function(obj) {
  if (is_atlas_runtime_bundle(obj)) {
    return(atlas_bundle_manifest(obj)[["features"]] %||% character(0))
  }
  rownames(obj)
}

atlas_cells <- function(obj) {
  if (is_atlas_runtime_bundle(obj)) {
    return(rownames(atlas_meta(obj)))
  }
  colnames(obj)
}

atlas_coords <- function(obj) {
  if (is_atlas_runtime_bundle(obj)) {
    return(atlas_bundle_umap_df(obj))
  }
  reduction_to_use <- preferred_umap_reduction(obj)
  validate(need(!is.null(reduction_to_use), "No compatible UMAP reduction was found in the selected Seurat object. Expected HarmonyUMAP or umap."))
  coords <- Embeddings(obj, reduction = reduction_to_use)
  validate(need(ncol(coords) >= 2, "The selected reduction does not have two dimensions."))
  coords <- as.data.frame(coords[, 1:2, drop = FALSE])
  colnames(coords) <- c("UMAP_1", "UMAP_2")
  coords
}

atlas_gene_manifest <- function(obj) {
  if (!is_atlas_runtime_bundle(obj)) {
    return(NULL)
  }
  manifest_path <- atlas_bundle_manifest(obj)[["gene_manifest_path"]] %||% NULL
  validate(need(!is.null(manifest_path) && file.exists(manifest_path), "Gene-expression manifest for the atlas runtime bundle was not found."))
  cached_read_rds(manifest_path)
}

atlas_gene_expression_path <- function(obj, gene) {
  validate(need(is_atlas_runtime_bundle(obj), "Per-gene atlas expression paths are only available for atlas runtime bundles."))
  manifest <- atlas_gene_manifest(obj)
  row <- manifest[manifest$gene %in% gene, , drop = FALSE]
  validate(need(nrow(row) == 1, paste0("Gene '", gene, "' is not available in the atlas runtime bundle.")))
  gene_dir <- atlas_bundle_manifest(obj)[["gene_dir"]] %||% NULL
  validate(need(!is.null(gene_dir) && dir.exists(gene_dir), "Gene-expression directory for the atlas runtime bundle was not found."))
  normalizePath(file.path(gene_dir, row$file_name[[1]]), winslash = "/", mustWork = FALSE)
}

atlas_expr_matrix <- function(obj) {
  if (is_atlas_runtime_bundle(obj)) {
    stop("Full expression matrices are not loaded at runtime for atlas bundles; use atlas_fetch_expression().", call. = FALSE)
  }
  get_expr_slot_safe(obj, atlas_default_assay(obj), "data")
}

atlas_fetch_expression <- function(obj, features, cells = NULL) {
  feature_ids <- unique(as.character(features))
  feature_ids <- feature_ids[feature_ids %in% atlas_features(obj)]
  validate(need(length(feature_ids) > 0, "None of the selected genes were found in the loaded atlas object."))

  if (!is_atlas_runtime_bundle(obj)) {
    expr <- get_expr_slot_safe(obj, atlas_default_assay(obj), "data")
    cell_ids <- cells %||% colnames(expr)
    cell_ids <- intersect(as.character(cell_ids), colnames(expr))
    validate(need(length(cell_ids) > 0, "No cells remain after filtering."))
    return(expr[feature_ids, cell_ids, drop = FALSE])
  }

  full_cells <- atlas_bundle_base_cells(obj)
  requested_cells <- cells %||% atlas_cells(obj)
  cell_ids <- intersect(as.character(requested_cells), full_cells)
  validate(need(length(cell_ids) > 0, "No cells remain after filtering."))

  use_all_cells <- identical(cell_ids, full_cells)
  cell_positions <- if (use_all_cells) seq_along(full_cells) else match(cell_ids, full_cells)
  position_map <- NULL
  if (!use_all_cells) {
    position_map <- integer(length(full_cells))
    position_map[cell_positions] <- seq_along(cell_positions)
  }

  out <- matrix(
    0,
    nrow = length(feature_ids),
    ncol = length(cell_ids),
    dimnames = list(feature_ids, cell_ids)
  )

  manifest <- atlas_gene_manifest(obj)
  manifest <- manifest[manifest$gene %in% feature_ids, , drop = FALSE]
  validate(need(nrow(manifest) > 0, "None of the selected genes were found in the atlas runtime bundle."))

  storage_mode <- atlas_bundle_manifest(obj)[["gene_storage"]] %||% "rds_sparse_block"
  if (identical(storage_mode, "rds_sparse_vector")) {
    for (i in seq_along(feature_ids)) {
      asset <- cached_read_rds(atlas_gene_expression_path(obj, feature_ids[[i]]))
      if (length(asset$index) == 0) {
        next
      }
      if (use_all_cells) {
        out[i, asset$index] <- asset$value
      } else {
        out_idx <- position_map[asset$index]
        keep <- out_idx > 0
        if (any(keep)) {
          out[i, out_idx[keep]] <- asset$value[keep]
        }
      }
    }
    return(out)
  }

  for (file_name in unique(manifest$file_name)) {
    block_path <- normalizePath(file.path(atlas_bundle_manifest(obj)[["gene_dir"]], file_name), winslash = "/", mustWork = FALSE)
    block_mat <- cached_read_rds(block_path)
    gene_rows <- manifest[manifest$file_name %in% file_name, , drop = FALSE]
    genes_in_block <- intersect(as.character(gene_rows$gene), rownames(block_mat))
    if (length(genes_in_block) == 0) {
      next
    }
    block_subset <- block_mat[genes_in_block, , drop = FALSE]
    if (!use_all_cells) {
      block_subset <- block_subset[, cell_positions, drop = FALSE]
    }
    row_idx <- match(genes_in_block, rownames(out))
    out[row_idx, seq_len(ncol(out))] <- as.matrix(block_subset)
  }

  out
}

atlas_subset_cells <- function(obj, cells) {
  cells <- intersect(as.character(cells), atlas_cells(obj))
  if (is_atlas_runtime_bundle(obj)) {
    manifest <- atlas_bundle_manifest(obj)
    manifest$meta_override <- atlas_meta(obj)[cells, , drop = FALSE]
    manifest$umap_override <- atlas_coords(obj)[cells, , drop = FALSE]
    manifest$n_cells <- nrow(manifest$meta_override)
    return(structure(manifest, class = "atlas_runtime_bundle"))
  }
  subset(obj, cells = cells)
}

atlas_set_meta_column <- function(obj, name, values) {
  if (is_atlas_runtime_bundle(obj)) {
    manifest <- atlas_bundle_manifest(obj)
    meta <- atlas_meta(obj)
    meta[[name]] <- values
    manifest$meta_override <- meta
    manifest$n_cells <- nrow(meta)
    return(structure(manifest, class = "atlas_runtime_bundle"))
  }
  obj[[name]] <- values
  obj
}

`[[.atlas_runtime_bundle` <- function(x, i, ..., drop = FALSE) {
  md <- atlas_meta(x)
  if (missing(i)) {
    return(md)
  }
  if (is.character(i) && length(i) == 1 && i %in% colnames(md)) {
    if (isTRUE(drop)) {
      return(md[[i]])
    }
    return(md[, i, drop = FALSE])
  }
  md[, i, drop = drop]
}

`$.atlas_runtime_bundle` <- function(x, name) {
  md <- atlas_meta(x)
  if (!is.null(md) && name %in% colnames(md)) {
    return(md[[name]])
  }
  atlas_bundle_manifest(x)[[name]]
}

seurat_cache_keys <- function() {
  ls(envir = .gc_seurat_object_cache, all.names = TRUE)
}

preload_seurat_object <- function(path) {
  invisible(safe_seurat_read(path))
}

cached_read_rds <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (exists(path, envir = .gc_rds_cache, inherits = FALSE)) {
    return(get(path, envir = .gc_rds_cache, inherits = FALSE))
  }
  value <- readRDS(path)
  assign(path, value, envir = .gc_rds_cache)
  value
}

prefer_app_slim_path <- function(path) {
  if (is.null(path) || !nzchar(path)) {
    return(path)
  }

  slim_path <- sub("\\.rds$", "_app_slim.rds", path, ignore.case = TRUE)
  if (!identical(slim_path, path) && file.exists(slim_path)) {
    return(slim_path)
  }
  path
}

matching_app_slim_path <- function(path) {
  if (is.null(path) || !nzchar(path)) {
    return(path)
  }
  sub("\\.rds$", "_app_slim.rds", path, ignore.case = TRUE)
}

matching_avg_cache_path <- function(path) {
  if (is.null(path) || !nzchar(path)) {
    return(path)
  }
  sub("\\.rds$", "_avg_cache.rds", path, ignore.case = TRUE)
}

matching_quadrant_cache_path <- function(path) {
  if (is.null(path) || !nzchar(path)) {
    return(path)
  }
  normalized_path <- normalized_source_path_for_mapping(path)
  if (!is.null(normalized_path) && identical(normalized_path, app_data_path("seurat_cancercells_final.rds"))) {
    merged_path <- prefer_app_slim_path(app_data_path("seurat_merged_TME_malignant_final_umap.rds"))
    return(sub("\\.rds$", "_quadrant_cache.rds", merged_path, ignore.case = TRUE))
  }
  sub("\\.rds$", "_quadrant_cache.rds", path, ignore.case = TRUE)
}

matching_quadrant_default_pdf_path <- function(path) {
  if (is.null(path) || !nzchar(path)) {
    return(path)
  }
  normalized_path <- normalized_source_path_for_mapping(path)
  if (!is.null(normalized_path) && identical(normalized_path, app_data_path("seurat_cancercells_final.rds"))) {
    merged_path <- prefer_app_slim_path(app_data_path("seurat_merged_TME_malignant_final_umap.rds"))
    return(sub("\\.rds$", "_quadrant_default.pdf", merged_path, ignore.case = TRUE))
  }
  sub("\\.rds$", "_quadrant_default.pdf", path, ignore.case = TRUE)
}

default_group_var <- function(obj) {
  preferred <- c("final_celltype", "celltype", "cell_type", "seurat_clusters", "ident", "orig.ident")
  available <- colnames(atlas_meta(obj))
  hit <- preferred[preferred %in% available]
  if (length(hit) > 0) {
    return(hit[[1]])
  }
  available[[1]] %||% NULL
}

available_reductions <- function(obj) {
  if (is_atlas_runtime_bundle(obj)) {
    return(atlas_reduction_name(obj))
  }
  names(obj@reductions)
}

preferred_umap_reduction <- function(obj) {
  if (is_atlas_runtime_bundle(obj)) {
    return(atlas_reduction_name(obj))
  }
  reductions <- available_reductions(obj)
  if ("HarmonyUMAP" %in% reductions) {
    return("HarmonyUMAP")
  }
  if ("umap" %in% reductions) {
    return("umap")
  }
  if ("UMAP" %in% reductions) {
    return("UMAP")
  }
  NULL
}

load_named_palette <- function(kind) {
  path <- color_rds_paths[[kind]]
  if (is.null(path) || !file.exists(path)) {
    return(setNames(character(0), character(0)))
  }

  pal <- readRDS(path)
  if (is.null(names(pal))) {
    return(setNames(as.character(pal), as.character(seq_along(pal))))
  }
  stats::setNames(as.character(pal), names(pal))
}

available_assays <- function(obj) {
  if (is_atlas_runtime_bundle(obj)) {
    return(atlas_default_assay(obj))
  }
  names(obj@assays)
}

available_features <- function(obj, source_path = NULL) {
  cache_key <- source_path %||% paste0("mem:", length(atlas_features(obj)), ":", nrow(atlas_meta(obj)))
  if (exists(cache_key, envir = .gc_feature_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .gc_feature_cache, inherits = FALSE))
  }
  features <- sort(unique(atlas_features(obj)))
  assign(cache_key, features, envir = .gc_feature_cache)
  features
}

object_summary <- function(obj, path = NULL, group_var = NULL) {
  meta <- atlas_meta(obj)
  active_group <- group_var %||% default_group_var(obj)
  active_levels <- if (!is.null(active_group) && active_group %in% colnames(meta)) {
    paste(unique(as.character(meta[[active_group]])), collapse = ", ")
  } else {
    ""
  }
  tibble::tibble(
    metric = c(
      "Source file",
      "Cells",
      "Features",
      "Assays",
      "Reductions",
      "Metadata columns",
      "Active identity levels"
    ),
    value = c(
      path %||% "Uploaded in session",
      format(nrow(meta), big.mark = ","),
      format(length(atlas_features(obj)), big.mark = ","),
      paste(available_assays(obj), collapse = ", "),
      paste(available_reductions(obj), collapse = ", "),
      paste(colnames(meta), collapse = ", "),
      active_levels
    )
  )
}

average_expression_long <- function(obj, genes, group_var, assay = NULL, slot = "data") {
  validate(need(length(genes) > 0, "Choose at least one gene."))
  validate(need(group_var %in% colnames(obj[[]]), "Selected grouping column is missing."))

  assay_to_use <- assay %||% DefaultAssay(obj)
  avg <- AverageExpression(
    object = obj,
    assays = assay_to_use,
    features = genes,
    group.by = group_var,
    slot = slot,
    verbose = FALSE
  )[[assay_to_use]]

  as.data.frame(as.table(as.matrix(avg)), stringsAsFactors = FALSE) |>
    dplyr::rename(gene = Var1, group = Var2, expression = Freq)
}

plot_average_heatmap <- function(avg_df) {
  ggplot(avg_df, aes(x = group, y = gene, fill = expression)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_viridis_c(option = "C") +
    labs(x = NULL, y = NULL, fill = "Average expr.") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
}

metadata_preview <- function(obj, n = 10) {
  head(atlas_meta(obj), n = n)
}

umap_field_labels <- c(
  final_celltype = "Cell Type",
  rev_pathological_subtype = "Pathological Group",
  rev_condition = "Primary/Normal",
  rev_molecular_subtype = "MSS/MSI",
  dataset = "Cohort",
  sampleID = "Sample ID",
  patientID = "Patient ID",
  TNM = "TNM",
  Phase = "Phase"
)

heatmap_filter_labels <- c(
  final_celltype = "Cell Type",
  dataset = "Cohort",
  rev_pathological_subtype = "Pathological Subtype",
  rev_condition = "Primary/Normal",
  rev_molecular_subtype = "MSS/MSI",
  sampleID = "Sample ID",
  patientID = "Patient ID",
  TNM = "TNM",
  Phase = "Phase"
)

umap_color_group_info <- function(obj, color_mode) {
  md <- atlas_meta(obj)

  validate(need(color_mode %in% names(umap_field_labels), "Unsupported UMAP color mode."))

  field_to_use <- color_mode
  if (identical(color_mode, "rev_pathological_subtype")) {
    validate(need("rev_pathological_subtype" %in% colnames(md), "Metadata column 'rev_pathological_subtype' was not found in the selected Seurat object."))
    field_to_use <- "rev_pathological_subtype"
  } else {
    validate(need(field_to_use %in% colnames(md), paste0("Metadata column '", field_to_use, "' was not found in the selected Seurat object.")))
  }

  values <- as.character(md[[field_to_use]])
  label <- unname(umap_field_labels[[color_mode]])

  if (identical(color_mode, "dataset")) {
    palette <- load_named_palette("dataset")
  } else if (identical(color_mode, "rev_pathological_subtype")) {
    palette <- load_named_palette("path_group")
  } else if (identical(color_mode, "final_celltype")) {
    if ("final_celltype" %in% colnames(md)) {
      values <- as.character(md$final_celltype)
      label <- unname(umap_field_labels[[color_mode]])
    } else {
      fallback_group <- default_group_var(obj)
      validate(need(!is.null(fallback_group) && fallback_group %in% colnames(md), "No cell-type-like grouping was found in the loaded atlas object."))
      values <- as.character(md[[fallback_group]])
      label <- "Cell Type (active.ident)"
    }
    palette <- load_named_palette("celltype")
  } else {
    palette <- setNames(character(0), character(0))
  }

  values[is.na(values) | !nzchar(values)] <- "NA"

  list(
    values = values,
    label = label,
    palette = palette
  )
}

palette_for_values <- function(values, named_palette) {
  values <- unique(as.character(values))
  values <- values[!is.na(values)]

  if (length(values) == 0) {
    return(setNames(character(0), character(0)))
  }

  out <- named_palette[intersect(values, names(named_palette))]
  missing_values <- setdiff(values, names(out))
  if (length(missing_values) > 0) {
    fallback <- grDevices::hcl.colors(length(missing_values), palette = "Dark 3")
    out <- c(out, stats::setNames(fallback, missing_values))
  }
  out
}

build_umap_plot_data <- function(obj, color_mode) {
  reduction_to_use <- preferred_umap_reduction(obj)
  validate(need(!is.null(reduction_to_use), "No compatible UMAP reduction was found in the loaded atlas object."))
  coords <- atlas_coords(obj)
  group_info <- umap_color_group_info(obj, color_mode)
  plot_df <- data.frame(
    UMAP_1 = coords$UMAP_1,
    UMAP_2 = coords$UMAP_2,
    group = group_info$values,
    stringsAsFactors = FALSE
  )

  palette <- palette_for_values(plot_df$group, group_info$palette)
  plot_df$group <- factor(plot_df$group, levels = names(palette))

  centers <- plot_df |>
    dplyr::group_by(group) |>
    dplyr::summarise(
      UMAP_1 = stats::median(UMAP_1, na.rm = TRUE),
      UMAP_2 = stats::median(UMAP_2, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    plot_df = plot_df,
    centers = centers,
    palette = palette,
    legend_title = group_info$label,
    reduction_name = reduction_to_use
  )
}

plot_harmony_umap <- function(obj, color_mode, label = TRUE, pt_size = 0.4) {
  plot_data <- build_umap_plot_data(obj, color_mode)

  p <- ggplot(plot_data$plot_df, aes(x = UMAP_1, y = UMAP_2, color = group)) +
    geom_point(size = pt_size, alpha = 0.8, shape = 16) +
    scale_color_manual(values = plot_data$palette, na.value = "#808080") +
    labs(
      x = paste0(plot_data$reduction_name, "_1"),
      y = paste0(plot_data$reduction_name, "_2"),
      color = plot_data$legend_title
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold")
    )

  if (isTRUE(label) && nrow(plot_data$centers) > 0) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = plot_data$centers,
        aes(x = UMAP_1, y = UMAP_2, label = group),
        inherit.aes = FALSE,
        size = 3,
        show.legend = FALSE
      )
    } else {
      p <- p + geom_text(
        data = plot_data$centers,
        aes(x = UMAP_1, y = UMAP_2, label = group),
        inherit.aes = FALSE,
        size = 3,
        show.legend = FALSE
      )
    }
  }

  p
}

xenium_color_labels <- c(
  celltype = "Cell Type",
  tumor_normal = "Tumor / Normal",
  neighborhood = "Neighborhood"
)

xenium_palette_for_var <- function(var_name, values) {
  values_chr <- unique(as.character(values))
  values_chr <- values_chr[!is.na(values_chr) & nzchar(values_chr)]

  if (length(values_chr) == 0) {
    return(setNames(character(0), character(0)))
  }

  if (identical(var_name, "celltype")) {
    return(palette_for_values(values_chr, load_named_palette("celltype")))
  }

  if (identical(var_name, "tumor_normal")) {
    base <- c(
      "Normal" = "#2CA25F",
      "Tumor" = "#D1495B"
    )
    return(palette_for_values(values_chr, base))
  }

  if (identical(var_name, "neighborhood")) {
    return(named_palette(values_chr, brewer_name = "Set2"))
  }

  named_palette(values_chr, brewer_name = "Set3")
}

xenium_image_is_plot_asset <- function(bundle) {
  source_image <- basename(bundle$image_meta$source_image %||% bundle$image_path %||% "")
  grepl("imagedimplot", source_image, ignore.case = TRUE)
}

xenium_sample_summary <- function(bundle) {
  tibble::tibble(
    metric = c(
      "Sample ID",
      "Cells",
      "Cell Types",
      "Neighborhoods",
      "Image Width",
      "Image Height"
    ),
    value = c(
      bundle$sample_id,
      format(nrow(bundle$cells), big.mark = ","),
      format(dplyr::n_distinct(bundle$cells$celltype), big.mark = ","),
      format(dplyr::n_distinct(bundle$cells$neighborhood), big.mark = ","),
      bundle$image_meta$displayed_width %||% bundle$image_meta$original_width %||% NA,
      bundle$image_meta$displayed_height %||% bundle$image_meta$original_height %||% NA
    )
  )
}

xenium_count_table <- function(cells_df, field) {
  validate(need(field %in% colnames(cells_df), paste0("'", field, "' is missing from the Xenium cell table.")))

  cells_df |>
    dplyr::count(.data[[field]], name = "n_cells", sort = TRUE) |>
    dplyr::mutate(
      fraction = n_cells / sum(n_cells),
      fraction = sprintf("%.1f%%", fraction * 100)
    ) |>
    dplyr::rename(category = 1)
}

xenium_overview_barplot <- function(cells_df, field, title) {
  count_df <- xenium_count_table(cells_df, field)
  count_df$category <- factor(count_df$category, levels = rev(count_df$category))

  ggplot(count_df, aes(x = category, y = n_cells, fill = category)) +
    geom_col(width = 0.8, show.legend = FALSE) +
    coord_flip() +
    scale_fill_manual(values = xenium_palette_for_var(field, as.character(count_df$category))) +
    labs(x = NULL, y = "Cells", title = title) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

plot_xenium_spatial_map <- function(bundle, color_by = "celltype", point_size = 0.45, alpha = 0.8, highlight_celltypes = character(0)) {
  validate(need(color_by %in% colnames(bundle$cells), paste0("'", color_by, "' is not available for Xenium coloring.")))

  plot_df <- bundle$cells |>
    dplyr::mutate(color_value = as.character(.data[[color_by]]))

  plot_df$color_value[is.na(plot_df$color_value) | !nzchar(plot_df$color_value)] <- "NA"
  highlight_celltypes <- unique(as.character(highlight_celltypes))
  highlight_celltypes <- highlight_celltypes[nzchar(highlight_celltypes) & !is.na(highlight_celltypes)]

  use_highlight_mode <- length(highlight_celltypes) > 0
  if (use_highlight_mode) {
    first_label <- highlight_celltypes[[1]]
    second_label <- if (length(highlight_celltypes) >= 2) highlight_celltypes[[2]] else NULL

    plot_df$plot_group <- "Other cells"
    plot_df$plot_group[plot_df$celltype %in% first_label] <- first_label
    if (!is.null(second_label) && !identical(second_label, first_label)) {
      plot_df$plot_group[plot_df$celltype %in% second_label] <- second_label
    }

    color_levels <- c("Other cells", first_label)
    color_values <- c(
      "Other cells" = "#C9CED6",
      stats::setNames("#2C7BE5", first_label)
    )

    if (!is.null(second_label) && !identical(second_label, first_label)) {
      color_levels <- c(color_levels, second_label)
      color_values <- c(color_values, stats::setNames("#D1495B", second_label))
    }

    plot_df$plot_group <- factor(plot_df$plot_group, levels = color_levels)
  } else {
    color_values <- xenium_palette_for_var(color_by, plot_df$color_value)
    plot_df$plot_group <- factor(plot_df$color_value, levels = names(color_values))
  }

  width <- as.numeric(bundle$image_meta$displayed_width %||% bundle$image_meta$original_width %||% max(plot_df$x_plot, na.rm = TRUE))
  height <- as.numeric(bundle$image_meta$displayed_height %||% bundle$image_meta$original_height %||% max(plot_df$y_plot, na.rm = TRUE))

  p <- ggplot()

  if (!xenium_image_is_plot_asset(bundle)) {
    validate(need(requireNamespace("png", quietly = TRUE), "Install the 'png' package to display Xenium background images."))
    img <- png::readPNG(bundle$image_path)
    p <- p + annotation_raster(img, xmin = 0, xmax = width, ymin = 0, ymax = height)
  }

  background_df <- plot_df[as.character(plot_df$plot_group) %in% "Other cells", , drop = FALSE]
  focus_df <- plot_df[!(as.character(plot_df$plot_group) %in% "Other cells"), , drop = FALSE]

  p +
    {
      if (nrow(background_df) > 0) {
        geom_point(
          data = background_df,
          aes(x = x_plot, y = y_plot, color = plot_group),
          size = point_size,
          alpha = if (use_highlight_mode) min(alpha, 0.25) else alpha,
          shape = 16
        )
      }
    } +
    {
      if (nrow(focus_df) > 0) {
        geom_point(
          data = focus_df,
          aes(x = x_plot, y = y_plot, color = plot_group),
          size = if (use_highlight_mode) point_size * 1.1 else point_size,
          alpha = alpha,
          shape = 16
        )
      }
    } +
    scale_color_manual(values = color_values, na.value = "#9AA5B1") +
    coord_fixed(xlim = c(0, width), ylim = c(0, height), expand = FALSE) +
    labs(color = if (use_highlight_mode) "Highlighted cell types" else (xenium_color_labels[[color_by]] %||% color_by)) +
    theme_void(base_size = 12) +
    theme(
      panel.background = element_rect(
        fill = if (xenium_image_is_plot_asset(bundle)) "#FBFCFD" else NA,
        color = NA
      ),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    ) +
    guides(color = guide_legend(override.aes = list(size = 5.5, alpha = 1)))
}

plot_xenium_neighborhood_abundance <- function(neighborhood_summary) {
  validate(need(nrow(neighborhood_summary) > 0, "Neighborhood summary is empty for this sample."))

  plot_df <- neighborhood_summary |>
    dplyr::arrange(dplyr::desc(fraction)) |>
    dplyr::mutate(neighborhood = factor(neighborhood, levels = rev(unique(neighborhood))))

  ggplot(plot_df, aes(x = neighborhood, y = fraction, fill = neighborhood)) +
    geom_col(width = 0.8, show.legend = FALSE) +
    coord_flip() +
    scale_fill_manual(values = xenium_palette_for_var("neighborhood", as.character(plot_df$neighborhood))) +
    scale_y_continuous(labels = function(x) sprintf("%d%%", round(x * 100))) +
    labs(x = NULL, y = "Fraction of cells", title = "Neighborhood abundance") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

plot_xenium_neighborhood_composition <- function(neighborhood_composition, neighborhood_summary = NULL) {
  validate(need(nrow(neighborhood_composition) > 0, "Neighborhood composition is empty for this sample."))

  plot_df <- neighborhood_composition

  neighborhood_levels <- if (!is.null(neighborhood_summary) && nrow(neighborhood_summary) > 0) {
    neighborhood_summary |>
      dplyr::arrange(dplyr::desc(fraction)) |>
      dplyr::pull(neighborhood) |>
      unique()
  } else {
    unique(plot_df$neighborhood)
  }

  celltype_levels <- plot_df |>
    dplyr::group_by(celltype) |>
    dplyr::summarise(total_fraction = sum(fraction, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(total_fraction)) |>
    dplyr::pull(celltype)

  plot_df$neighborhood <- factor(plot_df$neighborhood, levels = rev(neighborhood_levels))
  plot_df$celltype <- factor(plot_df$celltype, levels = celltype_levels)

  ggplot(plot_df, aes(x = celltype, y = neighborhood, fill = fraction)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_viridis_c(option = "C", labels = function(x) sprintf("%d%%", round(x * 100))) +
    labs(x = NULL, y = NULL, fill = "Fraction", title = "Neighborhood composition by cell type") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

plot_xenium_marker_overlay <- function(bundle, gene, point_size = 0.45, alpha = 0.8) {
  expr <- xenium_read_marker_expression(bundle, gene)
  plot_df <- bundle$cells
  plot_df$expr <- as.numeric(expr)
  plot_df$expr_plot <- log1p(plot_df$expr)

  width <- as.numeric(bundle$image_meta$displayed_width %||% bundle$image_meta$original_width %||% max(plot_df$x_plot, na.rm = TRUE))
  height <- as.numeric(bundle$image_meta$displayed_height %||% bundle$image_meta$original_height %||% max(plot_df$y_plot, na.rm = TRUE))

  p <- ggplot()

  if (!xenium_image_is_plot_asset(bundle)) {
    validate(need(requireNamespace("png", quietly = TRUE), "Install the 'png' package to display Xenium background images."))
    img <- png::readPNG(bundle$image_path)
    p <- p + annotation_raster(img, xmin = 0, xmax = width, ymin = 0, ymax = height)
  }

  zero_df <- plot_df[plot_df$expr_plot <= 0 | is.na(plot_df$expr_plot), , drop = FALSE]
  nonzero_df <- plot_df[plot_df$expr_plot > 0 & !is.na(plot_df$expr_plot), , drop = FALSE]
  if (nrow(nonzero_df) > 0) {
    nonzero_df <- nonzero_df[order(nonzero_df$expr_plot), , drop = FALSE]
  }

  p +
    {
      if (nrow(zero_df) > 0) {
        geom_point(
          data = zero_df,
          aes(x = x_plot, y = y_plot),
          color = "grey82",
          size = point_size,
          alpha = min(alpha, 0.35),
          shape = 16
        )
      }
    } +
    {
      if (nrow(nonzero_df) > 0) {
        geom_point(
          data = nonzero_df,
          aes(x = x_plot, y = y_plot, color = expr_plot),
          size = point_size,
          alpha = alpha,
          shape = 16
        )
      }
    } +
    scale_color_gradientn(colors = viridisLite::plasma(256), name = "log1p(count)") +
    coord_fixed(xlim = c(0, width), ylim = c(0, height), expand = FALSE) +
    labs(title = gene) +
    theme_void(base_size = 12) +
    theme(
      panel.background = element_rect(
        fill = if (xenium_image_is_plot_asset(bundle)) "#FBFCFD" else NA,
        color = NA
      ),
      plot.title = element_text(face = "bold", hjust = 0.02),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )
}

plot_feature_expression <- function(obj, features, pt_size = 0.35, xlim = NULL, ylim = NULL) {
  features <- unique(as.character(features))
  features <- features[features %in% atlas_features(obj)]
  validate(need(length(features) > 0, "None of the selected genes were found in the loaded atlas object."))

  coords <- atlas_coords(obj)
  md <- atlas_meta(obj)

  label_values <- if ("final_celltype" %in% colnames(md)) {
    as.character(md$final_celltype)
  } else {
    fallback_group <- default_group_var(obj)
    validate(need(!is.null(fallback_group) && fallback_group %in% colnames(md), "No grouping column is available for feature-plot labels."))
    as.character(md[[fallback_group]])
  }
  label_values[is.na(label_values) | !nzchar(label_values)] <- "NA"

  label_centers <- data.frame(
    UMAP_1 = coords$UMAP_1,
    UMAP_2 = coords$UMAP_2,
    label = label_values,
    stringsAsFactors = FALSE
  ) |>
    dplyr::group_by(label) |>
    dplyr::summarise(
      UMAP_1 = stats::median(UMAP_1, na.rm = TRUE),
      UMAP_2 = stats::median(UMAP_2, na.rm = TRUE),
      .groups = "drop"
    )

  if (!is.null(xlim) && length(xlim) == 2) {
    label_centers <- label_centers[
      is.na(label_centers$UMAP_1) | (label_centers$UMAP_1 >= min(xlim) & label_centers$UMAP_1 <= max(xlim)),
      ,
      drop = FALSE
    ]
  }
  if (!is.null(ylim) && length(ylim) == 2) {
    label_centers <- label_centers[
      is.na(label_centers$UMAP_2) | (label_centers$UMAP_2 >= min(ylim) & label_centers$UMAP_2 <= max(ylim)),
      ,
      drop = FALSE
    ]
  }

  expr <- as.data.frame(t(as.matrix(atlas_fetch_expression(obj, features = features, cells = rownames(md)))), stringsAsFactors = FALSE)

  plot_df <- do.call(
    rbind,
    lapply(features, function(feature) {
      values <- as.numeric(expr[[feature]])
      data.frame(
        UMAP_1 = coords$UMAP_1,
        UMAP_2 = coords$UMAP_2,
        feature = feature,
        expression = values,
        expression_nonzero = ifelse(values > 0, values, NA_real_),
        stringsAsFactors = FALSE
      )
    })
  )

  plot_df_zero <- plot_df[is.na(plot_df$expression_nonzero), , drop = FALSE]
  plot_df_nonzero <- plot_df[!is.na(plot_df$expression_nonzero), , drop = FALSE]
  if (nrow(plot_df_nonzero) > 0) {
    plot_df_nonzero <- plot_df_nonzero[order(plot_df_nonzero$expression_nonzero), , drop = FALSE]
  }

  p <- ggplot() +
    geom_point(
      data = plot_df_zero,
      aes(x = UMAP_1, y = UMAP_2),
      color = "grey80",
      size = pt_size,
      alpha = 0.7
    ) +
    geom_point(
      data = plot_df_nonzero,
      aes(x = UMAP_1, y = UMAP_2, color = expression_nonzero),
      size = pt_size,
      alpha = 0.9
    ) +
    scale_color_gradientn(
      colors = viridisLite::plasma(256),
      name = "Expression"
    ) +
    facet_wrap(~feature, ncol = min(2, length(features))) +
    labs(
      x = "UMAP_1",
      y = "UMAP_2"
    ) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold"),
      strip.text = element_blank(),
      strip.background = element_blank()
    )

  if (nrow(label_centers) > 0) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = label_centers,
        aes(x = UMAP_1, y = UMAP_2, label = label),
        inherit.aes = FALSE,
        size = 4.1,
        color = "black",
        box.padding = 0.25,
        segment.color = "grey50",
        show.legend = FALSE
      )
    } else {
      p <- p + geom_text(
        data = label_centers,
        aes(x = UMAP_1, y = UMAP_2, label = label),
        inherit.aes = FALSE,
        size = 4.1,
        color = "black",
        show.legend = FALSE
      )
    }
  }

  p
}

violin_display_label <- function(var_name) {
  labels <- c(
    ".violin_pathology_group" = "Pathological Group",
    ".violin_gc_status" = "Normal vs Matched GC",
    "rev_condition" = "Primary/Normal",
    "rev_pathological_subtype" = "Pathological Subtype",
    "rev_molecular_subtype" = "MSS/MSI",
    "dataset" = "Cohort",
    "final_celltype" = "Cell Type"
  )
  hit <- labels[var_name]
  if (length(hit) == 0 || is.na(hit)) {
    return(var_name)
  }
  unname(hit)
}

violin_factor_levels <- function(values, var_name) {
  values_chr <- as.character(values)
  observed <- unique(values_chr[!is.na(values_chr) & nzchar(values_chr)])

  preferred <- switch(
    var_name,
    ".violin_pathology_group" = c("Diffuse", "Intestinal", "Mixed"),
    ".violin_gc_status" = c("Normal", "Diffuse GC", "Intestinal GC", "Mixed GC"),
    "rev_condition" = c("Normal", "Primary", "GC"),
    "rev_molecular_subtype" = c("MSS", "MSI"),
    NULL
  )

  if (is.null(preferred)) {
    return(observed)
  }

  c(preferred[preferred %in% observed], sort(setdiff(observed, preferred)))
}

violin_palette_for_var <- function(var_name, values) {
  values_chr <- unique(as.character(values))
  values_chr <- values_chr[!is.na(values_chr) & nzchar(values_chr)]

  if (identical(var_name, ".violin_pathology_group") || identical(var_name, "rev_pathological_subtype")) {
    return(palette_for_values(values_chr, load_named_palette("path_group")))
  }
  if (identical(var_name, "dataset")) {
    return(palette_for_values(values_chr, load_named_palette("dataset")))
  }
  if (identical(var_name, "final_celltype")) {
    return(palette_for_values(values_chr, load_named_palette("celltype")))
  }
  if (identical(var_name, ".violin_gc_status") || identical(var_name, "rev_condition")) {
    path_pal <- load_named_palette("path_group")
    base <- c(
      "Normal" = unname(path_pal[["Normal"]] %||% "#4DAF4A"),
      "GC" = "#D95F02",
      "Primary" = "#D95F02",
      "Diffuse GC" = unname(path_pal[["Diffuse"]] %||% "#E41A1C"),
      "Intestinal GC" = unname(path_pal[["Intestinal"]] %||% "#377EB8"),
      "Mixed GC" = unname(path_pal[["Mixed"]] %||% "#984EA3")
    )
    return(palette_for_values(values_chr, base))
  }
  if (identical(var_name, "rev_molecular_subtype")) {
    base <- c(
      "MSS" = "#1B9E77",
      "MSI" = "#D95F02"
    )
    return(palette_for_values(values_chr, base))
  }

  named_palette(values_chr, brewer_name = "Set2")
}

GeomSplitViolin <- ggplot2::ggproto(
  "GeomSplitViolin",
  ggplot2::GeomViolin,
  draw_group = function(self, data, panel_params, coord, ...) {
    data <- transform(
      data,
      xminv = x - violinwidth * (x - xmin),
      xmaxv = x + violinwidth * (xmax - x)
    )

    is_left <- data[1, "group"] %% 2 == 1
    ordered <- data[order(data$y), , drop = FALSE]
    ordered$x <- if (is_left) ordered$xminv else ordered$xmaxv

    polygon <- rbind(
      data.frame(x = data$x[1], y = min(data$y), ordered[1, setdiff(names(ordered), c("x", "y")), drop = FALSE]),
      ordered,
      data.frame(x = data$x[1], y = max(data$y), ordered[nrow(ordered), setdiff(names(ordered), c("x", "y")), drop = FALSE]),
      data.frame(x = data$x[1], y = min(data$y), ordered[1, setdiff(names(ordered), c("x", "y")), drop = FALSE])
    )

    ggplot2::GeomPolygon$draw_panel(polygon, panel_params, coord, ...)
  }
)

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity",
                              ..., trim = TRUE, scale = "width", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomSplitViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(trim = trim, scale = scale, na.rm = na.rm, ...)
  )
}

plot_violin_expression <- function(obj, features, group_var, split_var = NULL) {
  features <- unique(as.character(features))
  features <- features[features %in% atlas_features(obj)]
  validate(need(length(features) > 0, "None of the selected genes were found in the loaded atlas object."))

  md <- atlas_meta(obj)
  validate(need(group_var %in% colnames(md), "Selected grouping column is not available."))
  if (!is.null(split_var)) {
    validate(need(split_var %in% colnames(md), "Selected split column is not available."))
  }
  expr <- as.matrix(atlas_fetch_expression(obj, features = features, cells = rownames(md)))
  expr <- t(expr)
  plot_data <- cbind(md[, unique(c(group_var, split_var)), drop = FALSE], as.data.frame(expr, stringsAsFactors = FALSE))
  plot_data$cell_id <- rownames(md)

  validate(need(group_var %in% colnames(plot_data), "Selected grouping column is not available."))
  plot_data[[group_var]] <- factor(
    as.character(plot_data[[group_var]]),
    levels = violin_factor_levels(plot_data[[group_var]], group_var)
  )

  if (!is.null(split_var)) {
    validate(need(split_var %in% colnames(plot_data), "Selected split column is not available."))
    plot_data[[split_var]] <- factor(
      as.character(plot_data[[split_var]]),
      levels = violin_factor_levels(plot_data[[split_var]], split_var)
    )
  }

  pieces <- lapply(features, function(feature) {
    df <- plot_data[, c("cell_id", group_var, if (!is.null(split_var)) split_var, feature), drop = FALSE]
    colnames(df)[ncol(df)] <- "expression"
    df$feature <- feature
    df
  })
  df_long <- do.call(rbind, pieces)
  df_long <- df_long[!is.na(df_long[[group_var]]) & nzchar(as.character(df_long[[group_var]])), , drop = FALSE]
  validate(need(nrow(df_long) > 0, "No cells remain for violin plotting after grouping/filtering."))

  if (identical(split_var, ".violin_gc_status")) {
    df_long <- df_long[!is.na(df_long[[split_var]]) & nzchar(as.character(df_long[[split_var]])), , drop = FALSE]
    validate(need(nrow(df_long) > 0, "No cells remain for violin plotting after split filtering."))

    df_long$.split_side <- ifelse(as.character(df_long[[split_var]]) == "Normal", "Normal", "GC")
    df_long$.split_side <- factor(df_long$.split_side, levels = c("Normal", "GC"))
    df_long$.split_group <- interaction(df_long[[group_var]], df_long$.split_side, lex.order = TRUE, drop = TRUE)

    fill_levels <- violin_factor_levels(df_long[[split_var]], split_var)
    fill_palette <- violin_palette_for_var(split_var, fill_levels)[fill_levels]

    return(
      ggplot(
        df_long,
        aes(
          x = .data[[group_var]],
          y = .data[["expression"]],
          fill = .data[[split_var]],
          group = .data[[".split_group"]]
        )
      ) +
        geom_split_violin(
          trim = TRUE,
          scale = "width",
          color = "grey20",
          linewidth = 0.2
        ) +
        scale_fill_manual(
          values = fill_palette,
          drop = FALSE,
          name = violin_display_label(split_var)
        ) +
        facet_wrap(~feature, ncol = min(2, length(features)), scales = "free_y") +
        labs(
          x = violin_display_label(group_var),
          y = "Expression Level"
        ) +
        theme_minimal(base_size = 12) +
        theme(
          panel.grid.minor = element_blank(),
          axis.title = element_text(face = "bold"),
          legend.position = "right"
        )
    )
  }

  if (is.null(split_var)) {
    group_levels <- violin_factor_levels(df_long[[group_var]], group_var)
    fill_palette <- violin_palette_for_var(group_var, group_levels)[group_levels]
    return(
      ggplot(
        df_long,
        aes(x = .data[[group_var]], y = .data[["expression"]], fill = .data[[group_var]])
      ) +
        geom_violin(trim = TRUE, scale = "width", color = "grey20", linewidth = 0.2) +
        scale_fill_manual(values = fill_palette, drop = FALSE, name = violin_display_label(group_var)) +
        facet_wrap(~feature, ncol = min(2, length(features)), scales = "free_y") +
        labs(
          x = violin_display_label(group_var),
          y = "Expression Level",
          fill = violin_display_label(group_var)
        ) +
        theme(
          axis.title = element_text(face = "bold"),
          legend.position = "right"
        )
    )
  }

  split_values <- as.character(df_long[[split_var]])
  split_levels <- violin_factor_levels(split_values, split_var)
  cols <- unname(violin_palette_for_var(split_var, split_levels)[split_levels])
  df_long[[split_var]] <- factor(as.character(df_long[[split_var]]), levels = split_levels)

  if (length(split_levels) <= 2) {
    df_long$.split_group <- interaction(df_long[[group_var]], df_long[[split_var]], lex.order = TRUE, drop = TRUE)
    return(
      ggplot(
        df_long,
        aes(
          x = .data[[group_var]],
          y = .data[["expression"]],
          fill = .data[[split_var]],
          group = .data[[".split_group"]]
        )
      ) +
        geom_split_violin(
          trim = TRUE,
          scale = "width",
          color = "grey20",
          linewidth = 0.2
        ) +
        scale_fill_manual(values = cols, drop = FALSE, name = violin_display_label(split_var)) +
        facet_wrap(~feature, ncol = min(2, length(features)), scales = "free_y") +
        labs(
          x = violin_display_label(group_var),
          y = "Expression Level",
          fill = violin_display_label(split_var)
        ) +
        theme_minimal(base_size = 12) +
        theme(
          axis.title = element_text(face = "bold"),
          legend.position = "right"
        )
    )
  }

  ggplot(
    df_long,
    aes(x = .data[[group_var]], y = .data[["expression"]], fill = .data[[split_var]])
  ) +
    geom_violin(trim = TRUE, scale = "width", color = "grey20", linewidth = 0.2, position = position_dodge(width = 0.85)) +
    scale_fill_manual(values = cols, drop = FALSE, name = violin_display_label(split_var)) +
    facet_wrap(~feature, ncol = min(2, length(features)), scales = "free_y") +
    labs(
      x = violin_display_label(group_var),
      y = "Expression Level",
      fill = violin_display_label(split_var)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

mapped_excel_path <- function(source_path) {
  mapping <- stats::setNames(
    c(
      app_data_path("markers_RNA_merged_AllCelltypes_final.xlsx"),
      app_data_path("TableS4_RNAmarkers_normalEpi.xlsx"),
      app_data_path("TableS5_RNAmarkers_Stromal.xlsx"),
      app_data_path("TableS6_RNAmarkers_myeloid.xlsx"),
      app_data_path("TableS7_RNAmarkers_CD8T.xlsx"),
      app_data_path("TableS8_RNAmarkers_CD4T.xlsx"),
      app_data_path("TableS9_RNAmarkers_B.xlsx"),
      app_data_path("MAST_DE_by_GCmucous_3malig.xlsx")
    ),
    c(
      app_data_path("seurat_merged_TME_malignant_final_umap.rds"),
      app_data_path("seurat_epithelial_normal_final_final.rds"),
      app_data_path("seurat_Stromal_final.rds"),
      app_data_path("seurat_Mye_final2.rds"),
      app_data_path("seurat_CD8T_final2.rds"),
      app_data_path("seurat_CD4T_final2.rds"),
      app_data_path("seurat_B_final2.rds"),
      app_data_path("seurat_cancercells_final.rds")
    )
  )

  if (is.null(source_path) || !nzchar(source_path)) {
    return(NULL)
  }

  normalized_path <- sub("_app_slim\\.rds$", ".rds", source_path, ignore.case = TRUE)
  if (!normalized_path %in% names(mapping)) {
    return(NULL)
  }
  mapping[[normalized_path]]
}

normalized_source_path_for_mapping <- function(source_path) {
  if (is.null(source_path) || !nzchar(source_path)) {
    return(NULL)
  }
  sub("_app_slim\\.rds$", ".rds", source_path, ignore.case = TRUE)
}

deg_display_mode <- function(source_path) {
  normalized_path <- normalized_source_path_for_mapping(source_path)
  if (is.null(normalized_path)) {
    return("quadrant")
  }

  special_paths <- c(
    app_data_path("seurat_merged_TME_malignant_final_umap.rds"),
    app_data_path("seurat_cancercells_final.rds")
  )

  if (normalized_path %in% special_paths) {
    return("quadrant")
  }

  "diffuse_intestinal"
}

excel_sheets_safe <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    return(character(0))
  }
  if (!requireNamespace("readxl", quietly = TRUE)) {
    return(character(0))
  }
  readxl::excel_sheets(path)
}

read_excel_sheet_safe <- function(path, sheet) {
  validate(need(!is.null(path) && file.exists(path), "Mapped Excel file was not found."))
  validate(need(requireNamespace("readxl", quietly = TRUE), "Install the 'readxl' package to show Excel results in the Overview tab."))
  readxl::read_excel(path, sheet = sheet)
}

derive_subtype4 <- function(final_group) {
  out <- rep(NA_character_, length(final_group))
  final_group <- as.character(final_group)

  out[grepl("_Normal$", final_group)] <- "Normal"
  out[grepl("^Diffuse", final_group)] <- "Diffuse"
  out[grepl("^Intestinal", final_group)] <- "Intestinal"
  out[grepl("^Mixed", final_group)] <- "Mixed"
  out
}

pick_sample_column <- function(md) {
  candidate_sample_cols <- c(
    "sample",
    "sample_id",
    "orig.ident",
    "patient_id",
    "patient",
    "Sample",
    "slide_id",
    "Slide",
    "dataset"
  )

  hits <- intersect(candidate_sample_cols, colnames(md))
  if (length(hits) == 0) {
    stop(
      "Could not find a sample column. Expected one of: ",
      paste(candidate_sample_cols, collapse = ", ")
    )
  }
  hits[[1]]
}

agg_sparse_mean <- function(X, groups_vec) {
  groups_vec <- as.character(groups_vec)
  cols_unique <- unique(groups_vec)
  j <- match(groups_vec, cols_unique)

  indicator <- Matrix::sparseMatrix(
    i = seq_along(groups_vec),
    j = j,
    x = 1,
    dims = c(length(groups_vec), length(cols_unique))
  )
  colnames(indicator) <- cols_unique

  sum_by_group <- X %*% indicator
  n_per_group <- Matrix::colSums(indicator)
  avg_mat <- sweep(sum_by_group, 2, n_per_group, "/")

  list(avg = avg_mat, groups = cols_unique, counts = n_per_group)
}

row_zscore <- function(M) {
  M <- as.matrix(M)
  storage.mode(M) <- "double"
  Z <- t(scale(t(M), center = TRUE, scale = TRUE))
  Z[!is.finite(Z)] <- 0
  Z
}

named_palette <- function(values, brewer_name = "Set3") {
  values <- unique(as.character(values))
  values <- values[!is.na(values)]

  if (length(values) == 0) {
    return(setNames(character(0), character(0)))
  }

  if (requireNamespace("RColorBrewer", quietly = TRUE) && length(values) <= 12) {
    max_n <- max(3, length(values))
    pal <- RColorBrewer::brewer.pal(max_n, brewer_name)[seq_along(values)]
  } else {
    pal <- grDevices::hcl.colors(length(values), palette = "Dynamic")
  }

  setNames(pal, values)
}

heatmap_gaps <- function(values) {
  values <- as.character(values)
  if (length(values) <= 1) {
    return(NULL)
  }

  runs <- rle(values)
  if (length(runs$lengths) <= 1) {
    return(NULL)
  }

  head(cumsum(runs$lengths), -1)
}

average_heatmap_cache_exists <- function(source_path) {
  cache_path <- matching_avg_cache_path(source_path)
  !is.null(cache_path) && file.exists(cache_path)
}

read_average_heatmap_cache <- function(source_path) {
  cache_path <- matching_avg_cache_path(source_path)
  validate(need(!is.null(cache_path) && file.exists(cache_path), "Average-expression cache was not found for this object. Build the cache first."))
  cache <- cached_read_rds(cache_path)
  validate(need(is.list(cache) && !is.null(cache$avg_mat) && !is.null(cache$meta), "Average-expression cache file is not in the expected format."))
  cache
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

available_heatmap_filter_fields <- function(cache) {
  fields <- names(heatmap_filter_labels)
  fields[fields %in% colnames(cache$meta)]
}

available_heatmap_display_fields <- function(cache) {
  c("sample_celltype", available_heatmap_filter_fields(cache))
}

heatmap_display_label <- function(field) {
  if (identical(field, "sample_celltype")) {
    return("Sample x Cell Type")
  }
  unname(heatmap_filter_labels[[field]] %||% field)
}

build_heatmap_display_choices <- function(cache) {
  fields <- available_heatmap_display_fields(cache)
  stats::setNames(fields, vapply(fields, heatmap_display_label, character(1)))
}

default_heatmap_display_field <- function(cache) {
  fields <- available_heatmap_display_fields(cache)
  preferred <- c("sample_celltype", "final_celltype", "rev_pathological_subtype", "dataset")
  hit <- preferred[preferred %in% fields]
  hit[[1]] %||% fields[[1]]
}

default_heatmap_genes <- function(features) {
  requested <- c("OLFM4", "ONCUT2", "CPS1", "CLDN18", "TFF1")
  resolved <- c()

  for (gene in requested) {
    if (gene %in% features) {
      resolved <- c(resolved, gene)
    } else if (identical(gene, "ONCUT2") && "ONECUT2" %in% features) {
      resolved <- c(resolved, "ONECUT2")
    }
  }

  hits <- unique(resolved)
  if (length(hits) > 0) {
    return(hits)
  }
  head(features, 4)
}

default_heatmap_filter_values <- function(cache, field) {
  if (is.null(cache) || is.null(cache$meta) || !field %in% colnames(cache$meta)) {
    return(character(0))
  }

  if (identical(field, "final_celltype")) {
    preferred <- c("malig_Diffuse", "malig_Mixed", "malig_Intestinal", "Metaplasia", "GC_mucous")
    values <- sort(unique(as.character(cache$meta[[field]])))
    return(preferred[preferred %in% values])
  }

  character(0)
}

filter_cached_heatmap_meta <- function(meta, filters) {
  keep <- rep(TRUE, nrow(meta))
  for (field in names(filters)) {
    values <- filters[[field]]
    if (length(values) > 0 && field %in% colnames(meta)) {
      keep <- keep & as.character(meta[[field]]) %in% as.character(values)
    }
  }
  keep
}

heatmap_annotation_colors <- function(annotation) {
  colors <- list()
  for (field in colnames(annotation)) {
    values <- as.character(annotation[[field]])
    if (field %in% c("final_celltype", "Cell Type")) {
      colors[[field]] <- palette_for_values(values, load_named_palette("celltype"))
    } else if (field %in% c("rev_pathological_subtype", "Pathological Subtype")) {
      colors[[field]] <- palette_for_values(values, load_named_palette("path_group"))
    } else if (field %in% c("dataset", "Cohort")) {
      colors[[field]] <- palette_for_values(values, load_named_palette("dataset"))
    } else if (field %in% c("rev_molecular_subtype", "MSS/MSI")) {
      colors[[field]] <- violin_palette_for_var("rev_molecular_subtype", values)
    } else {
      colors[[field]] <- named_palette(values, brewer_name = "Set3")
    }
  }
  colors
}

order_heatmap_columns <- function(meta) {
  preferred <- c(
    "final_celltype",
    "rev_pathological_subtype",
    "rev_molecular_subtype",
    "dataset",
    "final_group",
    "pseudobulk_id"
  )
  fields <- preferred[preferred %in% colnames(meta)]
  if (length(fields) == 0) {
    return(seq_len(nrow(meta)))
  }
  do.call(order, unname(lapply(fields, function(field) as.character(meta[[field]]))))
}

build_annotation_from_cache_meta <- function(meta) {
  annotation_fields <- c(
    "final_celltype",
    "rev_pathological_subtype",
    "rev_molecular_subtype",
    "dataset"
  )
  fields <- annotation_fields[annotation_fields %in% colnames(meta)]
  if (length(fields) == 0) {
    return(data.frame(row.names = rownames(meta)))
  }
  out <- meta[, fields, drop = FALSE]
  rownames(out) <- rownames(meta)
  colnames(out) <- vapply(fields, heatmap_display_label, character(1))
  out
}

aggregate_cached_matrix <- function(mat, group_values) {
  agg <- agg_sparse_mean(mat, group_values)
  mat_out <- as.matrix(agg$avg)
  colnames(mat_out) <- agg$groups
  mat_out
}

build_average_heatmap_data_live <- function(
  obj,
  genes,
  mode = c("celltype_group", "sample", "metadata"),
  selected_celltypes = NULL,
  group_var = NULL
) {
  mode <- match.arg(mode)
  validate(need(length(genes) > 0, "Choose at least one gene."))

  md <- atlas_meta(obj)

  if ("final_celltype" %in% colnames(md) && length(selected_celltypes %||% character(0)) > 0) {
    md <- md[md$final_celltype %in% selected_celltypes, , drop = FALSE]
  }

  validate(need(nrow(md) > 0, "No cells remain after filtering."))

  genes_found <- intersect(genes, atlas_features(obj))
  validate(need(length(genes_found) > 0, "None of the selected genes are present in the Seurat object."))
  expr <- atlas_fetch_expression(obj, features = genes_found, cells = rownames(md))

  subtype_colors <- c(
    Diffuse = "#E41A1C",
    Intestinal = "#377EB8",
    Mixed = "#984EA3",
    Normal = "#4DAF4A"
  )

  if (mode == "sample") {
    validate(need("final_group" %in% colnames(md), "The object must contain a 'final_group' column for sample heatmaps."))

    md$Subtype4 <- derive_subtype4(md$final_group)
    md <- md[!is.na(md$Subtype4), , drop = FALSE]
    expr <- expr[, rownames(md), drop = FALSE]
    validate(need(ncol(expr) > 0, "No cells remain after applying the final_group-based subtype filter."))

    sample_col <- pick_sample_column(md)
    agg <- agg_sparse_mean(expr, md[[sample_col]])
    mat <- agg$avg
    colnames(mat) <- agg$groups

    if ("dataset" %in% colnames(md)) {
      ann_samp <- md |>
        dplyr::mutate(.sample_key = .data[[sample_col]]) |>
        dplyr::group_by(.sample_key, final_group, Subtype4, dataset, .add = FALSE) |>
        dplyr::summarise(n = dplyr::n(), .groups = "drop_last")
    } else {
      ann_samp <- md |>
        dplyr::mutate(.sample_key = .data[[sample_col]], dataset = NA_character_) |>
        dplyr::group_by(.sample_key, final_group, Subtype4, dataset, .add = FALSE) |>
        dplyr::summarise(n = dplyr::n(), .groups = "drop_last")
    }

    ann_samp <- ann_samp |>
      dplyr::slice_max(n, n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::select(.sample_key, final_group, Subtype4, dataset) |>
      dplyr::distinct() |>
      as.data.frame()

    rownames(ann_samp) <- ann_samp$.sample_key
    ann_samp$.sample_key <- NULL
    ann_samp <- ann_samp[colnames(mat), , drop = FALSE]

    ord <- do.call(order, list(ann_samp$Subtype4, ann_samp$final_group, rownames(ann_samp)))
    mat <- mat[, ord, drop = FALSE]
    ann_samp <- ann_samp[ord, , drop = FALSE]

    final_group_colors <- named_palette(ann_samp$final_group, brewer_name = "Set3")
    ann_colors <- list(Subtype4 = subtype_colors, final_group = final_group_colors)
    if ("dataset" %in% colnames(ann_samp) && any(!is.na(ann_samp$dataset))) {
      ann_colors$dataset <- named_palette(ann_samp$dataset, brewer_name = "Set1")
    }

    annotation <- ann_samp[, intersect(c("Subtype4", "final_group", "dataset"), colnames(ann_samp)), drop = FALSE]
    title <- sprintf("Per-sample mean expression (genes=%d), z-scored by gene", nrow(mat))
    gaps <- heatmap_gaps(annotation$Subtype4)
  } else if (mode == "celltype_group") {
    validate(need("final_celltype" %in% colnames(md), "The object must contain a 'final_celltype' column for the Atlas heatmap."))
    validate(need("final_group" %in% colnames(md), "The object must contain a 'final_group' column for the Atlas heatmap."))

    md$Subtype4 <- derive_subtype4(md$final_group)
    md <- md[!is.na(md$Subtype4), , drop = FALSE]
    expr <- expr[, rownames(md), drop = FALSE]
    validate(need(ncol(expr) > 0, "No cells remain after applying the final_group-based subtype filter."))

    colkey <- interaction(md$final_celltype, md$final_group, drop = TRUE, sep = " | ")
    agg <- agg_sparse_mean(expr, colkey)
    mat <- agg$avg
    colnames(mat) <- agg$groups

    parts <- strsplit(colnames(mat), " | ", fixed = TRUE)
    ann_cf <- data.frame(
      final_celltype = vapply(parts, `[`, "", 1),
      final_group = vapply(parts, `[`, "", 2),
      stringsAsFactors = FALSE
    )
    rownames(ann_cf) <- colnames(mat)
    ann_cf$Subtype4 <- derive_subtype4(ann_cf$final_group)
    ann_cf <- ann_cf[colnames(mat), , drop = FALSE]

    ord <- do.call(order, list(ann_cf$Subtype4, ann_cf$final_celltype, ann_cf$final_group))
    mat <- mat[, ord, drop = FALSE]
    ann_cf <- ann_cf[ord, , drop = FALSE]

    final_group_colors <- named_palette(ann_cf$final_group, brewer_name = "Set3")
    annotation <- ann_cf[, c("Subtype4", "final_group"), drop = FALSE]
    ann_colors <- list(Subtype4 = subtype_colors, final_group = final_group_colors)
    title <- "Mean expression per (final_celltype x final_group), z-scored by gene"
    gaps <- heatmap_gaps(annotation$Subtype4)
  } else {
    validate(need(!is.null(group_var) && nzchar(group_var), "Choose a metadata column for the heatmap."))
    validate(need(group_var %in% colnames(md), "Selected metadata column is missing."))

    agg <- agg_sparse_mean(expr, md[[group_var]])
    mat <- agg$avg
    colnames(mat) <- agg$groups

    annotation <- data.frame(group = colnames(mat), stringsAsFactors = FALSE)
    rownames(annotation) <- annotation$group
    ann_colors <- list(group = named_palette(annotation$group, brewer_name = "Set3"))
    title <- sprintf("Mean expression per %s, z-scored by gene", group_var)
    gaps <- NULL
  }

  mat_z <- row_zscore(mat)
  avg_long <- as.data.frame(as.table(as.matrix(mat)), stringsAsFactors = FALSE)
  colnames(avg_long) <- c("gene", "group", "average_expression")

  list(
    matrix_z = mat_z,
    matrix_avg = as.matrix(mat),
    annotation = annotation,
    annotation_colors = ann_colors,
    gaps = gaps,
    title = title,
    long = avg_long
  )
}

build_average_heatmap_data <- function(
  genes,
  source_path = NULL,
  obj = NULL,
  display_by = "sample_celltype",
  filters = list()
) {
  validate(need(length(genes) > 0, "Choose at least one gene."))

  if (!is.null(source_path) && average_heatmap_cache_exists(source_path)) {
    cache <- read_average_heatmap_cache(source_path)
    meta <- cache$meta
    keep <- filter_cached_heatmap_meta(meta, filters)
    validate(need(any(keep), "No cached pseudobulk columns remain after filtering."))

    meta_sub <- meta[keep, , drop = FALSE]
    mat_sub <- cache$avg_mat[, rownames(meta_sub), drop = FALSE]

    genes_found <- intersect(genes, rownames(mat_sub))
    validate(need(length(genes_found) > 0, "None of the selected genes are present in the cached matrix."))
    mat_sub <- mat_sub[genes_found, , drop = FALSE]

    if (identical(display_by, "sample_celltype")) {
      ord <- order_heatmap_columns(meta_sub)
      meta_sub <- meta_sub[ord, , drop = FALSE]
      mat_out <- as.matrix(mat_sub[, rownames(meta_sub), drop = FALSE])
      annotation <- build_annotation_from_cache_meta(meta_sub)
      title <- sprintf("Cached mean expression per sample x cell type (genes=%d)", nrow(mat_out))
      gaps <- if ("final_celltype" %in% colnames(annotation)) heatmap_gaps(annotation$final_celltype) else NULL
      display_colnames <- FALSE
    } else {
      validate(need(display_by %in% colnames(meta_sub), "Selected display grouping is not available in the cache metadata."))
      group_values <- as.character(meta_sub[[display_by]])
      keep_group <- !is.na(group_values) & nzchar(group_values)
      validate(need(any(keep_group), "No cached columns have a value for the selected display grouping."))

      meta_sub <- meta_sub[keep_group, , drop = FALSE]
      mat_sub <- mat_sub[, rownames(meta_sub), drop = FALSE]
      group_values <- as.character(meta_sub[[display_by]])

      mat_out <- aggregate_cached_matrix(mat_sub, group_values)
      annotation <- data.frame(group = colnames(mat_out), stringsAsFactors = FALSE)
      rownames(annotation) <- annotation$group
      title <- sprintf("Cached mean expression grouped by %s (genes=%d)", heatmap_display_label(display_by), nrow(mat_out))
      gaps <- NULL
      display_colnames <- TRUE
    }

    annotation_colors <- heatmap_annotation_colors(annotation)
    mat_z <- row_zscore(mat_out)
    avg_long <- as.data.frame(as.table(as.matrix(mat_out)), stringsAsFactors = FALSE)
    colnames(avg_long) <- c("gene", "group", "average_expression")

    return(list(
      matrix_z = mat_z,
      matrix_avg = as.matrix(mat_out),
      annotation = annotation,
      annotation_colors = annotation_colors,
      gaps = gaps,
      title = title,
      long = avg_long,
      from_cache = TRUE,
      show_colnames = display_colnames
    ))
  }

  validate(need(!is.null(obj), "Average-expression cache was not found, and no Seurat object was available for fallback computation."))
  build_average_heatmap_data_live(
    obj = obj,
    genes = genes,
    mode = "celltype_group",
    selected_celltypes = filters$final_celltype %||% NULL,
    group_var = display_by
  )
}

plot_average_heatmap_from_script <- function(heatmap_data) {
  validate(need(requireNamespace("pheatmap", quietly = TRUE), "Install the 'pheatmap' package to use the average heatmap panel."))

  mat <- heatmap_data$matrix_z
  max_abs <- suppressWarnings(max(abs(range(mat, finite = TRUE))))
  if (!is.finite(max_abs) || max_abs == 0) {
    max_abs <- 1
  }

  pal_cont <- grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(200)
  p <- pheatmap::pheatmap(
    mat,
    color = pal_cont,
    breaks = seq(-max_abs, max_abs, length.out = 201),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_col = heatmap_data$annotation,
    annotation_colors = heatmap_data$annotation_colors,
    gaps_col = heatmap_data$gaps,
    show_colnames = isTRUE(heatmap_data$show_colnames %||% TRUE),
    angle_col = 90,
    fontsize_row = 12,
    silent = TRUE,
    main = heatmap_data$title
  )

  grid::grid.newpage()
  grid::grid.draw(p$gtable)
  invisible(p)
}

get_expr_slot_safe <- function(obj, assay_name, slot_name) {
  if (is_atlas_runtime_bundle(obj)) {
    stop("Full expression matrices are not available at runtime for atlas bundles.", call. = FALSE)
  }
  out <- tryCatch(
    GetAssayData(obj, assay = assay_name, layer = slot_name),
    error = function(e) NULL
  )
  if (!is.null(out)) {
    return(out)
  }
  GetAssayData(obj, assay = assay_name, slot = slot_name)
}

calc_sample_expression_pct <- function(obj, genes, assay_name, slot_name, threshold, sample_col, keep_cells = NULL) {
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (length(genes) == 0) {
    return(tibble::tibble(gene = character(), pct_samples_expressed = numeric()))
  }

  genes_present <- intersect(genes, atlas_features(obj))
  if (length(genes_present) == 0) {
    return(tibble::tibble(gene = genes, pct_samples_expressed = 0))
  }

  meta <- atlas_meta(obj)
  validate(need(sample_col %in% colnames(meta), paste0("Metadata column '", sample_col, "' is required for sample-percentage calculations.")))
  sample_ids <- meta[[sample_col]]
  names(sample_ids) <- rownames(meta)
  valid_cells <- rownames(meta)
  if (!is.null(keep_cells)) {
    valid_cells <- intersect(valid_cells, as.character(keep_cells))
  }
  sample_ids <- sample_ids[valid_cells]
  valid_keep <- !is.na(sample_ids) & nzchar(as.character(sample_ids))
  valid_cells <- valid_cells[valid_keep]
  sample_ids <- sample_ids[valid_keep]
  if (length(valid_cells) == 0) {
    return(tibble::tibble(gene = genes, pct_samples_expressed = 0))
  }

  expr <- atlas_fetch_expression(obj, features = genes_present, cells = valid_cells)
  sample_ids <- as.character(sample_ids)
  sample_levels <- sort(unique(sample_ids))
  sample_factor <- factor(sample_ids, levels = sample_levels)
  sample_index <- as.integer(sample_factor)
  n_cells_per_sample <- tabulate(sample_index, nbins = length(sample_levels))

  membership <- Matrix::sparseMatrix(
    i = seq_along(sample_index),
    j = sample_index,
    x = 1,
    dims = c(length(sample_index), length(sample_levels))
  )
  sample_sums <- expr %*% membership
  sample_means <- sweep(as.matrix(sample_sums), 2, n_cells_per_sample, "/")
  pct <- rowMeans(sample_means > threshold) * 100

  out <- tibble::tibble(gene = genes_present, pct_samples_expressed = as.numeric(pct))
  tibble::tibble(gene = genes) |>
    dplyr::left_join(out, by = "gene") |>
    dplyr::mutate(pct_samples_expressed = dplyr::coalesce(pct_samples_expressed, 0))
}

calc_quadrant_specific_sample_pct <- function(obj, fc_tbl, assay_name, slot_name, threshold, sample_col, group_col) {
  meta <- atlas_meta(obj)
  validate(need(group_col %in% colnames(meta), paste0("Metadata column '", group_col, "' is required for quadrant-specific sample percentages.")))

  group_values <- as.character(meta[[group_col]])
  cell_ids <- rownames(meta)

  cell_sets <- list(
    "Q1 (+/+) GC_all" = cell_ids[grepl("_Primary$", group_values)],
    "Q2 (-/+) Intestinal" = cell_ids[grepl("^Intestinal.*_Primary$", group_values)],
    "Q3 (-/-) Normal" = cell_ids[grepl("_Normal$", group_values)],
    "Q4 (+/-) Diffuse" = cell_ids[grepl("^Diffuse.*_Primary$", group_values)]
  )

  pct_tables <- lapply(names(cell_sets), function(quadrant_name) {
    pct_tbl <- calc_sample_expression_pct(
      obj = obj,
      genes = fc_tbl$gene,
      assay_name = assay_name,
      slot_name = slot_name,
      threshold = threshold,
      sample_col = sample_col,
      keep_cells = cell_sets[[quadrant_name]]
    )
    colnames(pct_tbl)[colnames(pct_tbl) == "pct_samples_expressed"] <- "pct_quadrant"
    pct_tbl$quadrant <- quadrant_name
    pct_tbl
  })

  dplyr::bind_rows(pct_tables) |>
    dplyr::right_join(dplyr::select(fc_tbl, gene, quadrant), by = c("gene", "quadrant")) |>
    dplyr::mutate(pct_samples_expressed = dplyr::coalesce(pct_quadrant, 0)) |>
    dplyr::select(gene, quadrant, pct_samples_expressed)
}

quadrant_sample_pct_from_cache <- function(sample_means_by_quadrant, fc_tbl, threshold) {
  validate(need(is.list(sample_means_by_quadrant), "Quadrant cache does not contain precomputed sample means."))

  pct_tables <- lapply(names(sample_means_by_quadrant), function(quadrant_name) {
    sample_means <- sample_means_by_quadrant[[quadrant_name]]
    genes <- fc_tbl$gene[fc_tbl$quadrant %in% quadrant_name]
    genes <- unique(as.character(genes))
    genes_present <- intersect(genes, rownames(sample_means))
    if (length(genes_present) == 0 || ncol(sample_means) == 0) {
      pct_tbl <- tibble::tibble(gene = genes, pct_quadrant = 0)
    } else {
      pct <- rowMeans(sample_means[genes_present, , drop = FALSE] > threshold) * 100
      pct_tbl <- tibble::tibble(gene = genes_present, pct_quadrant = as.numeric(pct))
      pct_tbl <- tibble::tibble(gene = genes) |>
        dplyr::left_join(pct_tbl, by = "gene") |>
        dplyr::mutate(pct_quadrant = dplyr::coalesce(pct_quadrant, 0))
    }
    pct_tbl$quadrant <- quadrant_name
    pct_tbl
  })

  dplyr::bind_rows(pct_tables) |>
    dplyr::right_join(dplyr::select(fc_tbl, gene, quadrant), by = c("gene", "quadrant")) |>
    dplyr::mutate(pct_samples_expressed = dplyr::coalesce(pct_quadrant, 0)) |>
    dplyr::select(gene, quadrant, pct_samples_expressed)
}

calc_group_fc <- function(
  obj,
  genes_to_use,
  assay_name,
  slot_name,
  celltype_col,
  group_col,
  path_group,
  primary_celltypes,
  normal_celltypes,
  pseudocount = 1e-3
) {
  meta <- atlas_meta(obj)
  sel <- meta[[group_col]] %in% paste0(path_group, c("_Primary", "_Normal"))
  if (!any(sel)) {
    return(NULL)
  }

  sub_obj <- atlas_subset_cells(obj, cells = rownames(meta)[sel])
  md <- atlas_meta(sub_obj)
  md$pn_group <- dplyr::case_when(
    md[[celltype_col]] %in% primary_celltypes ~ "Primary",
    md[[celltype_col]] %in% normal_celltypes ~ "Normal",
    TRUE ~ NA_character_
  )
  keep <- !is.na(md$pn_group)
  if (!any(keep)) {
    return(NULL)
  }

  md <- md[keep, , drop = FALSE]
  expr <- atlas_fetch_expression(sub_obj, features = genes_to_use, cells = rownames(md))

  primary_cells <- md$pn_group == "Primary"
  normal_cells <- md$pn_group == "Normal"
  if (sum(primary_cells) == 0 || sum(normal_cells) == 0) {
    return(NULL)
  }

  primary_mean <- Matrix::rowMeans(expr[, primary_cells, drop = FALSE])
  normal_mean <- Matrix::rowMeans(expr[, normal_cells, drop = FALSE])

  genes_present <- intersect(genes_to_use, rownames(expr))
  primary_mean <- primary_mean[genes_present]
  normal_mean <- normal_mean[genes_present]

  fc <- log2((primary_mean + pseudocount) / (normal_mean + pseudocount))
  tibble::tibble(
    gene = genes_present,
    !!paste0(tolower(path_group), "_fc") := fc
  )
}

build_quadrant_fc_cache <- function(
  obj,
  genes_to_use = NULL,
  assay_name = NULL,
  slot_name = "data",
  celltype_col = "final_celltype",
  group_col = "final_group",
  pseudocount = 1e-3
) {
  meta <- atlas_meta(obj)
  validate(need(celltype_col %in% colnames(meta), paste0("Metadata column '", celltype_col, "' is not available.")))
  validate(need(group_col %in% colnames(meta), paste0("Metadata column '", group_col, "' is not available.")))

  assay_to_use <- assay_name %||% atlas_default_assay(obj)
  genes_to_use <- unique(as.character(genes_to_use %||% atlas_features(obj)))
  genes_to_use <- genes_to_use[!is.na(genes_to_use) & nzchar(genes_to_use)]

  diffuse_fc <- calc_group_fc(
    obj = obj,
    genes_to_use = genes_to_use,
    assay_name = assay_to_use,
    slot_name = slot_name,
    celltype_col = celltype_col,
    group_col = group_col,
    path_group = "Diffuse",
    primary_celltypes = c("malig_Diffuse"),
    normal_celltypes = c("GC_mucous"),
    pseudocount = pseudocount
  )
  intestinal_fc <- calc_group_fc(
    obj = obj,
    genes_to_use = genes_to_use,
    assay_name = assay_to_use,
    slot_name = slot_name,
    celltype_col = celltype_col,
    group_col = group_col,
    path_group = "Intestinal",
    primary_celltypes = c("malig_Intestinal"),
    normal_celltypes = c("GC_mucous"),
    pseudocount = pseudocount
  )

  validate(need(!is.null(diffuse_fc) && !is.null(intestinal_fc), "Could not compute quadrant FC cache; required Diffuse/Intestinal malignant and GC_mucous cells were not found."))

  fc_tbl <- diffuse_fc |>
    dplyr::full_join(intestinal_fc, by = "gene") |>
    dplyr::filter(is.finite(diffuse_fc), is.finite(intestinal_fc)) |>
    dplyr::mutate(
      quadrant = dplyr::case_when(
        diffuse_fc >= 0 & intestinal_fc >= 0 ~ "Q1 (+/+) GC_all",
        diffuse_fc < 0 & intestinal_fc >= 0 ~ "Q2 (-/+) Intestinal",
        diffuse_fc < 0 & intestinal_fc < 0 ~ "Q3 (-/-) Normal",
        diffuse_fc >= 0 & intestinal_fc < 0 ~ "Q4 (+/-) Diffuse",
        TRUE ~ NA_character_
      ),
      distance = sqrt(diffuse_fc^2 + intestinal_fc^2)
    )

  sample_col <- default_quadrant_sample_col(obj)
  sample_means_by_quadrant <- if (!is.null(sample_col) && sample_col %in% colnames(meta)) {
    expr_all <- get_expr_slot_safe(obj, assay_to_use, slot_name)
    group_values <- as.character(meta[[group_col]])
    cell_ids <- rownames(meta)
    cell_sets <- list(
      "Q1 (+/+) GC_all" = cell_ids[grepl("_Primary$", group_values)],
      "Q2 (-/+) Intestinal" = cell_ids[grepl("^Intestinal.*_Primary$", group_values)],
      "Q3 (-/-) Normal" = cell_ids[grepl("_Normal$", group_values)],
      "Q4 (+/-) Diffuse" = cell_ids[grepl("^Diffuse.*_Primary$", group_values)]
    )

    lapply(cell_sets, function(keep_cells) {
      valid_cells <- intersect(keep_cells, rownames(meta))
      if (length(valid_cells) == 0) {
        return(matrix(numeric(0), nrow = length(genes_to_use), ncol = 0, dimnames = list(genes_to_use, character(0))))
      }
      sample_ids <- as.character(meta[valid_cells, sample_col, drop = TRUE])
      keep_sample <- !is.na(sample_ids) & nzchar(sample_ids)
      valid_cells <- valid_cells[keep_sample]
      sample_ids <- sample_ids[keep_sample]
      if (length(valid_cells) == 0) {
        return(matrix(numeric(0), nrow = length(genes_to_use), ncol = 0, dimnames = list(genes_to_use, character(0))))
      }

      expr_subset <- expr_all[genes_to_use, valid_cells, drop = FALSE]
      sample_levels <- sort(unique(sample_ids))
      sample_factor <- factor(sample_ids, levels = sample_levels)
      sample_index <- as.integer(sample_factor)
      n_cells_per_sample <- tabulate(sample_index, nbins = length(sample_levels))
      membership <- Matrix::sparseMatrix(
        i = seq_along(sample_index),
        j = sample_index,
        x = 1,
        dims = c(length(sample_index), length(sample_levels))
      )
      sample_sums <- expr_subset %*% membership
      sample_means <- sweep(as.matrix(sample_sums), 2, n_cells_per_sample, "/")
      rownames(sample_means) <- rownames(expr_subset)
      colnames(sample_means) <- sample_levels
      sample_means
    })
  } else {
    NULL
  }

  list(
    fc_tbl = fc_tbl,
    assay_name = assay_to_use,
    slot_name = slot_name,
    celltype_col = celltype_col,
    group_col = group_col,
    sample_col = sample_col,
    sample_means_by_quadrant = sample_means_by_quadrant,
    pseudocount = pseudocount,
    created_at = as.character(Sys.time())
  )
}

quadrant_cache_exists <- function(source_path) {
  cache_path <- matching_quadrant_cache_path(source_path)
  !is.null(cache_path) && file.exists(cache_path)
}

read_quadrant_cache <- function(source_path) {
  cache_path <- matching_quadrant_cache_path(source_path)
  validate(need(!is.null(cache_path) && file.exists(cache_path), "Quadrant scatter cache was not found for this object. Build the cache first."))
  cache <- cached_read_rds(cache_path)
  validate(need(is.list(cache) && !is.null(cache$fc_tbl), "Quadrant cache file is not in the expected format."))
  cache
}

matching_diffuse_intestinal_deg_cache_path <- function(source_path) {
  if (is.null(source_path) || !nzchar(source_path)) {
    return(NULL)
  }
  sub("\\.rds$", "_diffuse_intestinal_deg_cache.rds", source_path, ignore.case = TRUE)
}

diffuse_intestinal_deg_cache_exists <- function(source_path) {
  cache_path <- matching_diffuse_intestinal_deg_cache_path(source_path)
  !is.null(cache_path) && file.exists(cache_path)
}

read_diffuse_intestinal_deg_cache <- function(source_path) {
  cache_path <- matching_diffuse_intestinal_deg_cache_path(source_path)
  validate(need(!is.null(cache_path) && file.exists(cache_path), "Diffuse-vs-Intestinal DEG cache was not found for this object. Build the cache first."))
  cache <- cached_read_rds(cache_path)
  validate(need(is.list(cache) && !is.null(cache$deg_table), "Diffuse-vs-Intestinal DEG cache file is not in the expected format."))
  cache
}

deg_cluster_column <- function(obj) {
  if ("final_celltype" %in% colnames(atlas_meta(obj))) {
    return("final_celltype")
  }
  NULL
}

derive_diffuse_intestinal_group <- function(final_group) {
  values <- as.character(final_group)
  out <- rep(NA_character_, length(values))
  out[grepl("^Diffuse_", values)] <- "Diffuse"
  out[grepl("^Intestinal_", values)] <- "Intestinal"
  out
}

build_diffuse_intestinal_deg_cache <- function(
  avg_cache,
  group_col = "final_group",
  cluster_col = "final_celltype",
  min_samples_per_group = 3
) {
  validate(need(is.list(avg_cache) && !is.null(avg_cache$avg_mat) && !is.null(avg_cache$meta), "Average-expression cache is not in the expected format."))

  avg_mat <- as.matrix(avg_cache$avg_mat)
  meta <- as.data.frame(avg_cache$meta, stringsAsFactors = FALSE)

  validate(need(group_col %in% colnames(meta), paste0("Metadata column '", group_col, "' is required in the average-expression cache.")))
  validate(need(ncol(avg_mat) == nrow(meta), "Average-expression cache dimensions are inconsistent."))

  if (!is.null(cluster_col) && cluster_col %in% colnames(meta)) {
    cluster_values <- as.character(meta[[cluster_col]])
  } else {
    cluster_col <- NULL
    cluster_values <- rep("All cells", nrow(meta))
  }

  compare_group <- derive_diffuse_intestinal_group(meta[[group_col]])
  keep <- !is.na(compare_group) & !is.na(cluster_values) & nzchar(cluster_values)
  validate(need(any(keep), "No Diffuse/Intestinal pseudobulk columns were found for this object."))

  avg_mat <- avg_mat[, keep, drop = FALSE]
  meta <- meta[keep, , drop = FALSE]
  cluster_values <- cluster_values[keep]
  compare_group <- compare_group[keep]

  cluster_levels <- sort(unique(cluster_values))
  results <- lapply(cluster_levels, function(cluster_name) {
    idx <- which(cluster_values %in% cluster_name & compare_group %in% c("Diffuse", "Intestinal"))
    if (length(idx) == 0) {
      return(NULL)
    }

    cluster_groups <- compare_group[idx]
    n_diffuse <- sum(cluster_groups %in% "Diffuse")
    n_intestinal <- sum(cluster_groups %in% "Intestinal")
    if (n_diffuse < min_samples_per_group || n_intestinal < min_samples_per_group) {
      return(NULL)
    }

    expr <- avg_mat[, idx, drop = FALSE]
    diffuse_idx <- which(cluster_groups %in% "Diffuse")
    intestinal_idx <- which(cluster_groups %in% "Intestinal")

    avg_log2fc <- rowMeans(expr[, diffuse_idx, drop = FALSE]) - rowMeans(expr[, intestinal_idx, drop = FALSE])
    pct_1 <- rowMeans(expr[, diffuse_idx, drop = FALSE] > 0)
    pct_2 <- rowMeans(expr[, intestinal_idx, drop = FALSE] > 0)

    p_val <- vapply(seq_len(nrow(expr)), function(i) {
      stats::wilcox.test(
        x = expr[i, diffuse_idx],
        y = expr[i, intestinal_idx],
        exact = FALSE
      )$p.value
    }, numeric(1))

    data.frame(
      gene = rownames(expr),
      avg_log2FC = as.numeric(avg_log2fc),
      pct.1 = as.numeric(pct_1),
      pct.2 = as.numeric(pct_2),
      p_val = as.numeric(p_val),
      p_val_adj = stats::p.adjust(p_val, method = "BH"),
      cluster = cluster_name,
      n_diffuse = n_diffuse,
      n_intestinal = n_intestinal,
      comparison = "Diffuse vs Intestinal",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  deg_table <- dplyr::bind_rows(results)
  validate(need(nrow(deg_table) > 0, "No valid Diffuse-vs-Intestinal DEG results were computed for this object."))

  list(
    deg_table = deg_table,
    group_col = group_col,
    cluster_col = cluster_col %||% "All cells",
    compare_groups = c("Diffuse", "Intestinal"),
    source_path = avg_cache$source_path %||% NULL,
    sample_col = avg_cache$sample_col %||% NULL,
    analysis_level = "sample_pseudobulk_wilcox",
    built_at = as.character(Sys.time())
  )
}

plot_diffuse_intestinal_volcano <- function(
  df,
  max_p_adj = 0.05,
  min_abs_log2fc = 0.25,
  top_n_labels = 15,
  max_abs_log2fc_display = 100
) {
  validate(need(nrow(df) > 0, "No DEG rows are available for the selected cell type."))
  validate(need("avg_log2FC" %in% colnames(df), "The DEG table is missing avg_log2FC."))
  validate(need("p_val_adj" %in% colnames(df), "The DEG table is missing p_val_adj."))

  pvals <- as.numeric(df$p_val_adj)
  smallest_nonzero <- suppressWarnings(min(pvals[pvals > 0], na.rm = TRUE))
  if (!is.finite(smallest_nonzero)) {
    smallest_nonzero <- 1e-300
  }

  plot_df <- df |>
    dplyr::filter(is.na(avg_log2FC) | abs(avg_log2FC) <= max_abs_log2fc_display) |>
    dplyr::mutate(
      p_val_adj_plot = dplyr::if_else(is.na(p_val_adj) | p_val_adj <= 0, smallest_nonzero, p_val_adj),
      neg_log10_p = -log10(p_val_adj_plot),
      significance = dplyr::case_when(
        !is.na(p_val_adj) & p_val_adj <= max_p_adj & avg_log2FC >= min_abs_log2fc ~ "Diffuse higher",
        !is.na(p_val_adj) & p_val_adj <= max_p_adj & avg_log2FC <= -min_abs_log2fc ~ "Intestinal higher",
        TRUE ~ "NS"
      )
    )

  validate(need(nrow(plot_df) > 0, paste0("No volcano points remain after excluding |log2FC| > ", max_abs_log2fc_display, ".")))

  label_df <- plot_df |>
    dplyr::filter(significance != "NS") |>
    dplyr::arrange(p_val_adj, dplyr::desc(abs(avg_log2FC))) |>
    dplyr::slice_head(n = top_n_labels)

  colors <- c(
    "Diffuse higher" = "#D1495B",
    "Intestinal higher" = "#2C7BE5",
    "NS" = "#B8C0CC"
  )

  p <- ggplot(plot_df, aes(x = avg_log2FC, y = neg_log10_p)) +
    geom_vline(xintercept = c(-min_abs_log2fc, min_abs_log2fc), linetype = "dashed", color = "grey70") +
    geom_hline(yintercept = -log10(max_p_adj), linetype = "dashed", color = "grey70") +
    geom_point(aes(color = significance), alpha = 0.8, size = 1.4) +
    scale_color_manual(values = colors, drop = FALSE) +
    labs(
      x = "Diffuse vs Intestinal (log2FC)",
      y = expression(-log[10]("adjusted p-value")),
      color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      legend.position = "top"
    )

  if (nrow(label_df) > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = label_df,
      aes(label = gene),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.35,
      show.legend = FALSE
    )
  }

  p
}

default_quadrant_sample_col <- function(obj) {
  if ("sampleID" %in% colnames(atlas_meta(obj))) {
    return("sampleID")
  }
  NULL
}

quadrant_sample_denominators <- function(obj, sample_col = "sampleID", group_col = "final_group") {
  meta <- atlas_meta(obj)
  validate(need(sample_col %in% colnames(meta), paste0("Metadata column '", sample_col, "' is required for the quadrant scatter.")))
  validate(need(group_col %in% colnames(meta), paste0("Metadata column '", group_col, "' is required for the quadrant scatter.")))

  group_values <- as.character(meta[[group_col]])
  sample_values <- as.character(meta[[sample_col]])

  count_samples <- function(mask) {
    vals <- unique(sample_values[mask & !is.na(sample_values) & nzchar(sample_values)])
    length(vals)
  }

  tibble::tibble(
    quadrant = c("Q1 (+/+) GC_all", "Q2 (-/+) Intestinal", "Q3 (-/-) Normal", "Q4 (+/-) Diffuse"),
    n_samples = c(
      count_samples(grepl("_Primary$", group_values)),
      count_samples(grepl("^Intestinal.*_Primary$", group_values)),
      count_samples(grepl("_Normal$", group_values)),
      count_samples(grepl("^Diffuse.*_Primary$", group_values))
    )
  )
}

build_quadrant_deg_data_from_cache <- function(
  obj,
  cache,
  sample_col = NULL,
  sample_expr_slot = NULL,
  sample_expr_threshold = 0.1,
  label_min_pct_samples = 25,
  label_qvalue = 0.05,
  force_label_genes = ""
) {
  base_data <- build_quadrant_deg_base_from_cache(
    obj = obj,
    cache = cache,
    sample_col = sample_col,
    sample_expr_slot = sample_expr_slot,
    sample_expr_threshold = sample_expr_threshold
  )

  apply_quadrant_label_settings(
    quadrant_base_data = base_data,
    label_min_pct_samples = label_min_pct_samples,
    label_qvalue = label_qvalue,
    force_label_genes = force_label_genes
  )
}

build_quadrant_deg_base_from_cache <- function(
  obj,
  cache,
  sample_col = NULL,
  sample_expr_slot = NULL,
  sample_expr_threshold = 0.1
) {
  validate(need(!is.null(cache$fc_tbl), "Quadrant cache is missing the FC table."))

  sample_col <- sample_col %||% default_quadrant_sample_col(obj)
  validate(need(!is.null(sample_col) && sample_col %in% colnames(atlas_meta(obj)), "A patient/sample column is required for the quadrant scatter."))

  sample_expr_slot <- sample_expr_slot %||% cache$slot_name %||% "data"

  plot_df <- cache$fc_tbl |>
    dplyr::filter(is.finite(diffuse_fc), is.finite(intestinal_fc))

  sample_pct_tbl <- if (!is.null(cache$sample_means_by_quadrant)) {
    quadrant_sample_pct_from_cache(
      sample_means_by_quadrant = cache$sample_means_by_quadrant,
      fc_tbl = plot_df,
      threshold = sample_expr_threshold
    )
  } else {
    assay_to_use <- cache$assay_name %||% atlas_default_assay(obj)
    calc_quadrant_specific_sample_pct(
      obj = obj,
      fc_tbl = plot_df,
      assay_name = assay_to_use,
      slot_name = sample_expr_slot,
      threshold = sample_expr_threshold,
      sample_col = sample_col,
      group_col = cache$group_col %||% "final_group"
    )
  }

  plot_df <- plot_df |>
    dplyr::left_join(sample_pct_tbl, by = c("gene", "quadrant")) |>
    dplyr::mutate(pct_samples_expressed = dplyr::coalesce(pct_samples_expressed, 0))

  mahalanobis_pool <- plot_df |>
    dplyr::filter(pct_samples_expressed > 0)

  if (nrow(mahalanobis_pool) >= 3) {
    coords <- as.matrix(mahalanobis_pool[, c("diffuse_fc", "intestinal_fc")])
    cov_mat <- stats::cov(coords)
    if (det(cov_mat) <= 0 || any(!is.finite(cov_mat))) {
      cov_mat <- stats::cov(coords + matrix(stats::rnorm(length(coords), sd = 1e-6), ncol = 2))
    }
    mdist <- stats::mahalanobis(coords, colMeans(coords), cov_mat)
    mahalanobis_tbl <- mahalanobis_pool |>
      dplyr::mutate(
        mahalanobis = mdist,
        p_value = 1 - stats::pchisq(mahalanobis, df = 2)
      ) |>
      dplyr::mutate(q_value_bh = stats::p.adjust(p_value, method = "BH")) |>
      dplyr::select(gene, mahalanobis, p_value, q_value_bh)
    plot_df <- plot_df |>
      dplyr::left_join(mahalanobis_tbl, by = "gene")
  } else {
    plot_df$mahalanobis <- NA_real_
    plot_df$p_value <- NA_real_
    plot_df$q_value_bh <- NA_real_
  }

  if (nrow(plot_df) > 1) {
    spearman <- suppressWarnings(stats::cor.test(plot_df$diffuse_fc, plot_df$intestinal_fc, method = "spearman"))
  } else {
    spearman <- NULL
  }

  list(
    plot_df = plot_df,
    spearman = spearman,
    size_legend_label = sprintf("%% relevant samples (mean %s > %.2g)", sample_expr_slot, sample_expr_threshold)
  )
}

apply_quadrant_label_settings <- function(
  quadrant_base_data,
  label_min_pct_samples = 25,
  label_qvalue = 0.05,
  force_label_genes = ""
) {
  plot_df <- quadrant_base_data$plot_df

  force_label_genes <- unlist(strsplit(force_label_genes %||% "", "[,;]"))
  force_label_genes <- trimws(force_label_genes)
  force_label_genes <- unique(force_label_genes[nzchar(force_label_genes)])

  label_df <- plot_df |>
    dplyr::filter(
      is.finite(q_value_bh),
      q_value_bh < label_qvalue,
      pct_samples_expressed >= label_min_pct_samples
    )

  if (length(force_label_genes) > 0) {
    label_df <- dplyr::bind_rows(
      label_df,
      dplyr::filter(plot_df, gene %in% force_label_genes)
    ) |>
      dplyr::distinct(gene, .keep_all = TRUE)
  }

  if (nrow(plot_df) > 1) {
    spearman <- quadrant_base_data$spearman
  } else {
    spearman <- NULL
  }

  sig_label_name <- sprintf("Significant (q<%s)", format(label_qvalue, digits = 2))
  plot_df <- plot_df |>
    dplyr::mutate(
      force_labeled = gene %in% force_label_genes,
      sig_label = dplyr::case_when(
        is.finite(q_value_bh) & q_value_bh < label_qvalue ~ sig_label_name,
        TRUE ~ "NS"
      ),
      point_color_group = dplyr::case_when(
        force_labeled ~ "Force-labeled",
        is.finite(q_value_bh) & q_value_bh < label_qvalue ~ sig_label_name,
        TRUE ~ "NS"
      )
    )

  list(
    plot_df = plot_df,
    label_df = label_df,
    spearman = spearman,
    sig_label_name = sig_label_name,
    size_legend_label = quadrant_base_data$size_legend_label
  )
}

build_quadrant_deg_data <- function(
  obj,
  genes_to_use,
  assay_name = NULL,
  slot_name = "data",
  sample_expr_slot = NULL,
  celltype_col = "final_celltype",
  group_col = "final_group",
  sample_col = "orig.ident",
  sample_expr_threshold = 0.1,
  label_min_pct_samples = 25,
  label_qvalue = 0.05,
  force_label_genes = "",
  pseudocount = 1e-3
) {
  validate(need(length(genes_to_use) > 0, "Choose at least one gene for the quadrant DEG plot."))
  validate(need(celltype_col %in% colnames(obj@meta.data), paste0("Metadata column '", celltype_col, "' is not available.")))
  validate(need(group_col %in% colnames(obj@meta.data), paste0("Metadata column '", group_col, "' is not available.")))
  validate(need(sample_col %in% colnames(obj@meta.data), paste0("Metadata column '", sample_col, "' is not available.")))

  assay_to_use <- assay_name %||% DefaultAssay(obj)
  sample_expr_slot <- sample_expr_slot %||% slot_name
  genes_to_use <- unique(as.character(genes_to_use))
  genes_to_use <- genes_to_use[!is.na(genes_to_use) & nzchar(genes_to_use)]

  diffuse_fc <- calc_group_fc(
    obj = obj,
    genes_to_use = genes_to_use,
    assay_name = assay_to_use,
    slot_name = slot_name,
    celltype_col = celltype_col,
    group_col = group_col,
    path_group = "Diffuse",
    primary_celltypes = c("malig_Diffuse"),
    normal_celltypes = c("GC_mucous"),
    pseudocount = pseudocount
  )
  intestinal_fc <- calc_group_fc(
    obj = obj,
    genes_to_use = genes_to_use,
    assay_name = assay_to_use,
    slot_name = slot_name,
    celltype_col = celltype_col,
    group_col = group_col,
    path_group = "Intestinal",
    primary_celltypes = c("malig_Intestinal"),
    normal_celltypes = c("GC_mucous"),
    pseudocount = pseudocount
  )

  validate(need(!is.null(diffuse_fc) && !is.null(intestinal_fc), "Quadrant DEG requires malignant-style Atlas labels including Diffuse/Intestinal final_group values and GC_mucous, malig_Diffuse, malig_Intestinal cell types."))

  fc_tbl <- diffuse_fc |>
    dplyr::full_join(intestinal_fc, by = "gene") |>
    dplyr::mutate(
      quadrant = dplyr::case_when(
        diffuse_fc >= 0 & intestinal_fc >= 0 ~ "Q1 (+/+)",
        diffuse_fc < 0 & intestinal_fc >= 0 ~ "Q2 (-/+)",
        diffuse_fc < 0 & intestinal_fc < 0 ~ "Q3 (-/-)",
        diffuse_fc >= 0 & intestinal_fc < 0 ~ "Q4 (+/-)",
        TRUE ~ NA_character_
      ),
      distance = sqrt(diffuse_fc^2 + intestinal_fc^2)
    )

  sample_pct_tbl <- calc_quadrant_specific_sample_pct(
    obj = obj,
    fc_tbl = fc_tbl,
    assay_name = assay_to_use,
    slot_name = sample_expr_slot,
    threshold = sample_expr_threshold,
    sample_col = sample_col,
    group_col = group_col
  )

  plot_df <- fc_tbl |>
    dplyr::filter(is.finite(diffuse_fc), is.finite(intestinal_fc)) |>
    dplyr::left_join(sample_pct_tbl, by = c("gene", "quadrant")) |>
    dplyr::mutate(pct_samples_expressed = dplyr::coalesce(pct_samples_expressed, 0))

  mahalanobis_pool <- plot_df |>
    dplyr::filter(pct_samples_expressed > 0)

  if (nrow(mahalanobis_pool) >= 3) {
    coords <- as.matrix(mahalanobis_pool[, c("diffuse_fc", "intestinal_fc")])
    cov_mat <- stats::cov(coords)
    if (det(cov_mat) <= 0 || any(!is.finite(cov_mat))) {
      cov_mat <- stats::cov(coords + matrix(stats::rnorm(length(coords), sd = 1e-6), ncol = 2))
    }
    mdist <- stats::mahalanobis(coords, colMeans(coords), cov_mat)
    mahalanobis_tbl <- mahalanobis_pool |>
      dplyr::mutate(
        mahalanobis = mdist,
        p_value = 1 - stats::pchisq(mahalanobis, df = 2)
      ) |>
      dplyr::mutate(q_value_bh = stats::p.adjust(p_value, method = "BH")) |>
      dplyr::select(gene, mahalanobis, p_value, q_value_bh)
    plot_df <- plot_df |>
      dplyr::left_join(mahalanobis_tbl, by = "gene")
  } else {
    plot_df$mahalanobis <- NA_real_
    plot_df$p_value <- NA_real_
    plot_df$q_value_bh <- NA_real_
  }

  force_label_genes <- unlist(strsplit(force_label_genes %||% "", "[,;]"))
  force_label_genes <- trimws(force_label_genes)
  force_label_genes <- unique(force_label_genes[nzchar(force_label_genes)])

  label_df <- plot_df |>
    dplyr::filter(
      is.finite(q_value_bh),
      q_value_bh < label_qvalue,
      pct_samples_expressed >= label_min_pct_samples
    )

  if (length(force_label_genes) > 0) {
    label_df <- dplyr::bind_rows(
      label_df,
      dplyr::filter(plot_df, gene %in% force_label_genes)
    ) |>
      dplyr::distinct(gene, .keep_all = TRUE)
  }

  if (nrow(plot_df) > 1) {
    spearman <- suppressWarnings(stats::cor.test(plot_df$diffuse_fc, plot_df$intestinal_fc, method = "spearman"))
  } else {
    spearman <- NULL
  }

  sig_label_name <- sprintf("Significant (q<%s)", format(label_qvalue, digits = 2))
  plot_df <- plot_df |>
    dplyr::mutate(
      force_labeled = gene %in% force_label_genes,
      sig_label = dplyr::case_when(
        is.finite(q_value_bh) & q_value_bh < label_qvalue ~ sig_label_name,
        TRUE ~ "NS"
      ),
      point_color_group = dplyr::case_when(
        force_labeled ~ "Force-labeled",
        is.finite(q_value_bh) & q_value_bh < label_qvalue ~ sig_label_name,
        TRUE ~ "NS"
      )
    )

  list(
    plot_df = plot_df,
    label_df = label_df,
    spearman = spearman,
    sig_label_name = sig_label_name,
    size_legend_label = sprintf("%% relevant samples (mean %s > %.2g)", sample_expr_slot, sample_expr_threshold)
  )
}

plot_quadrant_deg <- function(quadrant_data) {
  validate(need(nrow(quadrant_data$plot_df) > 0, "No finite fold-change pairs were available for the selected genes."))

  plot_df <- quadrant_data$plot_df
  label_df <- quadrant_data$label_df
  force_label_genes <- plot_df |>
    dplyr::filter(.data$force_labeled) |>
    dplyr::pull(.data$gene)
  force_point_df <- plot_df |>
    dplyr::filter(.data$force_labeled)
  force_label_df <- label_df |>
    dplyr::filter(.data$gene %in% force_label_genes)
  regular_label_df <- label_df |>
    dplyr::filter(!(.data$gene %in% force_label_df$gene))
  sig_label_name <- quadrant_data$sig_label_name
  size_legend_label <- quadrant_data$size_legend_label
  spearman <- quadrant_data$spearman
  color_levels <- c("Force-labeled", sig_label_name, "NS")
  color_values <- c(
    "Force-labeled" = "#1B9E77",
    stats::setNames("red", sig_label_name),
    "NS" = "#2c7fb8"
  )
  x_range <- range(plot_df$diffuse_fc, finite = TRUE)
  y_range <- range(plot_df$intestinal_fc, finite = TRUE)
  x_span <- diff(x_range)
  y_span <- diff(y_range)
  if (!is.finite(x_span) || x_span == 0) x_span <- 1
  if (!is.finite(y_span) || y_span == 0) y_span <- 1
  quadrant_labels <- tibble::tibble(
    x = c(
      x_range[2] - 0.03 * x_span,
      x_range[1] + 0.03 * x_span,
      x_range[1] + 0.03 * x_span,
      x_range[2] - 0.03 * x_span
    ),
    y = c(
      y_range[2] - 0.015 * y_span,
      y_range[2] - 0.03 * y_span,
      y_range[1] + 0.03 * y_span,
      y_range[1] + 0.03 * y_span
    ),
    label = c(
      "Global GC (Q1)",
      "Intestinal GC (Q2)",
      "Global Normal (Q3)",
      "Diffuse GC (Q4)"
    ),
    hjust = c(1, 0, 0, 1),
    vjust = c(1, 1, 0, 0)
  )

  p <- ggplot(plot_df, aes(x = diffuse_fc, y = intestinal_fc)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_text(
      data = quadrant_labels,
      aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
      inherit.aes = FALSE,
      size = 5.2,
      fontface = "bold",
      color = "black",
      alpha = 0.9
    ) +
    geom_point(aes(color = point_color_group, size = pct_samples_expressed), alpha = 0.7) +
    geom_point(
      data = force_point_df,
      aes(size = pct_samples_expressed),
      shape = 21,
      stroke = 1.2,
      fill = NA,
      color = "#1B9E77",
      show.legend = FALSE
    ) +
    labs(
      x = "Diffuse (Normal vs Primary, log2FC)",
      y = "Intestinal (Normal vs Primary, log2FC)",
      color = NULL,
      size = size_legend_label
    ) +
    scale_color_manual(values = color_values, breaks = color_levels, drop = FALSE) +
    scale_size_area(max_size = 5, limits = c(0, 100)) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      legend.position = "top"
    )

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    if (nrow(regular_label_df) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = regular_label_df,
        aes(label = gene),
        size = 2.8,
        max.overlaps = Inf,
        box.padding = 0.3,
        show.legend = FALSE
      )
    }
    if (nrow(force_label_df) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = force_label_df,
        aes(label = gene),
        color = "#1B9E77",
        size = 4.8,
        fontface = "bold",
        max.overlaps = Inf,
        box.padding = 0.45,
        point.padding = 0.3,
        segment.color = "#1B9E77",
        min.segment.length = 0,
        show.legend = FALSE
      )
    }
  }

  if (!is.null(spearman)) {
    p_label <- if (is.finite(spearman$p.value) && spearman$p.value < 0.001) {
      "P-value < 0.001"
    } else {
      sprintf("P-value = %.2g", spearman$p.value)
    }
    p <- p +
      annotate(
        "text",
        x = -Inf,
        y = Inf,
        label = sprintf("Spearman rho = %.2f", spearman$estimate),
        hjust = -0.1,
        vjust = 1.1,
        size = 3.2
      ) +
      annotate(
        "text",
        x = -Inf,
        y = Inf,
        label = p_label,
        hjust = -0.1,
        vjust = 2.4,
        size = 3.2
      )
  }

  p
}
