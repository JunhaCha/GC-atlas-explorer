#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) {
    return(args[[idx + 1]])
  }
  default
}

input_path <- get_arg("--input")
output_path <- get_arg("--output")
cache_path <- get_arg("--cache")
pdf_width <- as.numeric(get_arg("--width", "10.5"))
pdf_height <- as.numeric(get_arg("--height", "7"))

if (is.null(input_path) || !nzchar(input_path)) {
  stop("Provide --input <seurat_app_slim.rds>")
}

if (is.null(output_path) || !nzchar(output_path)) {
  output_path <- sub("\\.rds$", "_quadrant_default.pdf", input_path, ignore.case = TRUE)
}

if (is.null(cache_path) || !nzchar(cache_path)) {
  cache_path <- sub("\\.rds$", "_quadrant_cache.rds", input_path, ignore.case = TRUE)
}

get_expr_slot_safe <- function(obj, assay_name = NULL, slot_name = "data") {
  assay_name <- assay_name %||% obj@active.assay
  assay_obj <- obj@assays[[assay_name]]
  if (inherits(assay_obj, "Assay5")) {
    return(SeuratObject::LayerData(obj, assay = assay_name, layer = slot_name))
  }
  slot(assay_obj, slot_name)
}

calc_sample_expression_pct <- function(obj, genes, assay_name, slot_name, threshold, sample_col, keep_cells = NULL) {
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (length(genes) == 0) {
    return(tibble::tibble(gene = character(), pct_samples_expressed = numeric()))
  }

  expr_all <- get_expr_slot_safe(obj, assay_name, slot_name)
  genes_present <- intersect(genes, rownames(expr_all))
  if (length(genes_present) == 0) {
    return(tibble::tibble(gene = genes, pct_samples_expressed = 0))
  }

  sample_ids <- obj@meta.data[[sample_col]]
  names(sample_ids) <- rownames(obj@meta.data)
  sample_ids <- sample_ids[colnames(expr_all)]
  valid_cells <- !is.na(sample_ids) & nzchar(as.character(sample_ids))
  if (!is.null(keep_cells)) {
    valid_cells <- valid_cells & colnames(expr_all) %in% keep_cells
  }
  if (!any(valid_cells)) {
    return(tibble::tibble(gene = genes, pct_samples_expressed = 0))
  }

  expr <- expr_all[genes_present, valid_cells, drop = FALSE]
  sample_ids <- as.character(sample_ids[valid_cells])
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
  group_values <- as.character(obj@meta.data[[group_col]])
  cell_ids <- rownames(obj@meta.data)

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

build_quadrant_deg_data_from_cache <- function(
  obj,
  cache,
  sample_col = "sampleID",
  sample_expr_slot = NULL,
  sample_expr_threshold = 0.1,
  label_min_pct_samples = 10,
  label_qvalue = 0.05
) {
  sample_expr_slot <- sample_expr_slot %||% cache$slot_name %||% "data"
  assay_to_use <- cache$assay_name %||% obj@active.assay

  plot_df <- cache$fc_tbl |>
    dplyr::filter(is.finite(diffuse_fc), is.finite(intestinal_fc))

  sample_pct_tbl <- calc_quadrant_specific_sample_pct(
    obj = obj,
    fc_tbl = plot_df,
    assay_name = assay_to_use,
    slot_name = sample_expr_slot,
    threshold = sample_expr_threshold,
    sample_col = sample_col,
    group_col = cache$group_col %||% "final_group"
  )

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
        p_value = 1 - stats::pchisq(mahalanobis, df = 2),
        q_value_bh = stats::p.adjust(p_value, method = "BH")
      ) |>
      dplyr::select(gene, mahalanobis, p_value, q_value_bh)
    plot_df <- plot_df |>
      dplyr::left_join(mahalanobis_tbl, by = "gene")
  } else {
    plot_df$mahalanobis <- NA_real_
    plot_df$p_value <- NA_real_
    plot_df$q_value_bh <- NA_real_
  }

  sig_label_name <- sprintf("Significant (q<%s)", format(label_qvalue, digits = 2))
  label_df <- plot_df |>
    dplyr::filter(
      is.finite(q_value_bh),
      q_value_bh < label_qvalue,
      pct_samples_expressed >= label_min_pct_samples
    )

  if (nrow(plot_df) > 1) {
    spearman <- suppressWarnings(stats::cor.test(plot_df$diffuse_fc, plot_df$intestinal_fc, method = "spearman"))
  } else {
    spearman <- NULL
  }

  plot_df <- plot_df |>
    dplyr::mutate(
      force_labeled = FALSE,
      point_color_group = dplyr::case_when(
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
  plot_df <- quadrant_data$plot_df
  label_df <- quadrant_data$label_df
  sig_label_name <- quadrant_data$sig_label_name
  size_legend_label <- quadrant_data$size_legend_label
  spearman <- quadrant_data$spearman

  p <- ggplot(plot_df, aes(x = diffuse_fc, y = intestinal_fc)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(aes(color = point_color_group, size = pct_samples_expressed), alpha = 0.7) +
    labs(
      x = "Diffuse (Normal vs Primary, log2FC)",
      y = "Intestinal (Normal vs Primary, log2FC)",
      color = NULL,
      size = size_legend_label
    ) +
    scale_color_manual(
      values = c(stats::setNames("red", sig_label_name), "NS" = "#2c7fb8"),
      breaks = c(sig_label_name, "NS"),
      drop = FALSE
    ) +
    scale_size_area(max_size = 5, limits = c(0, 100)) +
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
      size = 2.8,
      max.overlaps = Inf,
      box.padding = 0.3,
      show.legend = FALSE
    )
  }

  if (!is.null(spearman)) {
    p <- p + annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = sprintf("Spearman rho = %.2f\np = %.2g", spearman$estimate, spearman$p.value),
      hjust = -0.1,
      vjust = 1.1,
      size = 3.2
    )
  }

  p
}

obj <- readRDS(input_path)
cache <- readRDS(cache_path)
quadrant_data <- build_quadrant_deg_data_from_cache(
  obj = obj,
  cache = cache,
  sample_col = "sampleID",
  sample_expr_slot = cache$slot_name %||% "data",
  sample_expr_threshold = 0.1,
  label_min_pct_samples = 10,
  label_qvalue = 0.05
)

p <- plot_quadrant_deg(quadrant_data)
grDevices::pdf(output_path, width = pdf_width, height = pdf_height)
print(p)
grDevices::dev.off()

cat("Wrote default quadrant PDF to ", output_path, "\n", sep = "")
