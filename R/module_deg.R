deg_ui <- function(id) {
  ns <- NS(id)

  fluidPage(
    fluidRow(
      column(
        width = 3,
        h5("Precomputed DEG Filters"),
        selectizeInput(
          ns("cluster_filter"),
          "Cluster filter",
          choices = NULL,
          multiple = TRUE,
          options = list(placeholder = "Leave empty to keep all clusters")
        ),
        selectizeInput(
          ns("gene_filter"),
          "Gene filter",
          choices = NULL,
          multiple = TRUE,
          options = list(placeholder = "Leave empty to keep all genes", maxOptions = 100000)
        ),
        numericInput(ns("max_p_adj"), "Max adjusted p-value", value = 0.05, min = 0, max = 1, step = 0.01),
        numericInput(ns("min_abs_log2fc"), "Min absolute log2FC", value = 0.25, min = 0, step = 0.05),
        numericInput(ns("min_pct_1"), "Min pct.1", value = 0.10, min = 0, max = 1, step = 0.05),
        numericInput(ns("min_pct_2"), "Min pct.2", value = 0, min = 0, max = 1, step = 0.05),
        selectInput(
          ns("sort_by"),
          "Sort table by",
          choices = c(
            "Adjusted p-value" = "p_val_adj",
            "Absolute log2FC" = "abs_log2FC",
            "Cluster" = "cluster",
            "Gene" = "gene"
          ),
          selected = "abs_log2FC"
        ),
        checkboxInput(ns("sort_desc"), "Sort descending", TRUE),
        br(),
        downloadButton(ns("download_deg"), "Download filtered DEG CSV"),
        tags$hr(),
        h5("Quadrant Scatter"),
        numericInput(ns("sample_expr_threshold"), "Sample expression threshold", value = 0.1, min = 0, step = 0.01),
        numericInput(ns("label_min_pct_samples"), "Min % samples for labels", value = 10, min = 0, max = 100, step = 5),
        numericInput(ns("label_qvalue"), "Mahalanobis q-value (BH)", value = 0.05, min = 0, max = 1, step = 0.005),
        textInput(ns("force_label_genes"), "Force-label genes", placeholder = "Comma-separated genes")
      ),
      column(
        width = 9,
        card(
          card_header("Wilcoxon DEG results"),
          textOutput(ns("deg_status")),
          br(),
          dataTableOutput(ns("deg_table"))
        ),
        br(),
        card(
          card_header("Normal epi vs malignant DEG (Intestinal and Diffuse)"),
          textOutput(ns("quadrant_status")),
          br(),
          plotOutput(ns("quadrant_plot"), height = 520)
        ),
        br(),
        card(
          card_header("Quadrant Results"),
          dataTableOutput(ns("quadrant_table"))
        )
      )
    )
  )
}

deg_sheet_choices <- function(path) {
  validate(need(!is.null(path) && file.exists(path), "No precomputed DEG workbook is mapped for this Seurat object yet."))
  validate(need(requireNamespace("readxl", quietly = TRUE), "Install the 'readxl' package to browse precomputed DEG workbooks."))

  sheets <- readxl::excel_sheets(path)
  validate(need(length(sheets) > 0, "The mapped DEG workbook does not contain any sheets."))
  c("All sheets" = "__all__", stats::setNames(sheets, sheets))
}

read_deg_sheet <- function(path, sheet) {
  readxl::read_excel(path, sheet = sheet) |>
    as.data.frame(stringsAsFactors = FALSE) |>
    dplyr::mutate(sheet = sheet, .before = 1)
}

read_precomputed_deg_table <- function(path, sheet_choice) {
  validate(need(!is.null(path) && file.exists(path), "No precomputed DEG workbook is mapped for this Seurat object yet."))
  validate(need(requireNamespace("readxl", quietly = TRUE), "Install the 'readxl' package to browse precomputed DEG workbooks."))

  sheet_map <- deg_sheet_choices(path)
  sheet_values <- unname(sheet_map)

  if (identical(sheet_choice, "__all__")) {
    sheets <- setdiff(sheet_values, "__all__")
    out <- lapply(sheets, function(sheet) read_deg_sheet(path, sheet))
    return(dplyr::bind_rows(out))
  }

  validate(need(sheet_choice %in% sheet_values, "Choose a valid workbook sheet."))
  read_deg_sheet(path, sheet_choice)
}

filter_precomputed_deg_table <- function(
  df,
  cluster_filter = character(0),
  gene_filter = character(0),
  max_p_adj = 0.05,
  min_abs_log2fc = 0,
  min_pct_1 = 0,
  min_pct_2 = 0,
  sort_by = "p_val_adj",
  sort_desc = FALSE
) {
  out <- df

  if (length(cluster_filter) > 0) {
    cluster_filter <- as.character(cluster_filter)
    if ("sheet" %in% colnames(out)) {
      out <- out[as.character(out$sheet) %in% cluster_filter, , drop = FALSE]
    } else if ("cluster" %in% colnames(out)) {
      out <- out[as.character(out$cluster) %in% cluster_filter, , drop = FALSE]
    }
  }
  if ("gene" %in% colnames(out) && length(gene_filter) > 0) {
    out <- out[as.character(out$gene) %in% as.character(gene_filter), , drop = FALSE]
  }
  if ("p_val_adj" %in% colnames(out) && !is.null(max_p_adj) && is.finite(max_p_adj)) {
    out <- out[is.na(out$p_val_adj) | out$p_val_adj <= max_p_adj, , drop = FALSE]
  }
  if ("avg_log2FC" %in% colnames(out) && !is.null(min_abs_log2fc) && is.finite(min_abs_log2fc)) {
    out <- out[is.na(out$avg_log2FC) | abs(out$avg_log2FC) >= min_abs_log2fc, , drop = FALSE]
  }
  if ("pct.1" %in% colnames(out) && !is.null(min_pct_1) && is.finite(min_pct_1)) {
    out <- out[is.na(out$`pct.1`) | out$`pct.1` >= min_pct_1, , drop = FALSE]
  }
  if ("pct.2" %in% colnames(out) && !is.null(min_pct_2) && is.finite(min_pct_2)) {
    out <- out[is.na(out$`pct.2`) | out$`pct.2` >= min_pct_2, , drop = FALSE]
  }

  if (identical(sort_by, "abs_log2FC") && "avg_log2FC" %in% colnames(out)) {
    out$abs_log2FC <- abs(out$avg_log2FC)
    out <- out[order(out$abs_log2FC, decreasing = isTRUE(sort_desc), na.last = TRUE), , drop = FALSE]
  } else if (sort_by %in% colnames(out)) {
    out <- out[order(out[[sort_by]], decreasing = isTRUE(sort_desc), na.last = TRUE), , drop = FALSE]
  }

  out
}

format_deg_table_for_display <- function(df) {
  if (nrow(df) == 0) {
    return(df)
  }

  out <- df

  if ("gene" %in% colnames(out) && !"genes" %in% colnames(out)) {
    out$genes <- out$gene
  }

  drop_cols <- intersect(c("sheet", "p_val", "abs_log2FC", "gene"), colnames(out))
  if (length(drop_cols) > 0) {
    out <- out[, setdiff(colnames(out), drop_cols), drop = FALSE]
  }

  preferred_order <- c("genes", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster")
  present_first <- preferred_order[preferred_order %in% colnames(out)]
  remaining <- setdiff(colnames(out), present_first)
  out[, c(present_first, remaining), drop = FALSE]
}

deg_server <- function(id, loaded) {
  moduleServer(id, function(input, output, session) {
    deg_workbook_path <- reactive({
      req(loaded()$source_path)
      mapped_excel_path(loaded()$source_path)
    })

    deg_table_raw <- reactive({
      path <- deg_workbook_path()
      validate(need(!is.null(path) && file.exists(path), "No precomputed DEG workbook is mapped for this Seurat object yet."))
      validate(need(requireNamespace("readxl", quietly = TRUE), "Install the 'readxl' package to browse precomputed DEG workbooks."))
      read_precomputed_deg_table(path, "__all__")
    })

    quadrant_cache <- reactive({
      req(loaded()$source_path)
      if (!quadrant_cache_exists(loaded()$source_path)) {
        return(NULL)
      }
      read_quadrant_cache(loaded()$source_path)
    })

    quadrant_base_results <- reactive({
      req(loaded()$obj)
      cache <- quadrant_cache()
      validate(need(!is.null(cache), "Quadrant scatter cache is not available for this object yet."))

      build_quadrant_deg_base_from_cache(
        obj = loaded()$obj,
        cache = cache,
        sample_expr_slot = cache$slot_name %||% "data",
        sample_expr_threshold = input$sample_expr_threshold
      )
    })

    quadrant_results <- reactive({
      req(quadrant_base_results())

      apply_quadrant_label_settings(
        quadrant_base_data = quadrant_base_results(),
        label_min_pct_samples = input$label_min_pct_samples,
        label_qvalue = input$label_qvalue,
        force_label_genes = input$force_label_genes
      )
    })

    observeEvent(deg_table_raw(), {
      df <- deg_table_raw()
      cluster_choices <- if ("sheet" %in% colnames(df)) {
        sort(unique(as.character(df$sheet)))
      } else if ("cluster" %in% colnames(df)) {
        sort(unique(as.character(df$cluster)))
      } else {
        character(0)
      }
      gene_choices <- if ("gene" %in% colnames(df)) sort(unique(as.character(df$gene))) else character(0)

      updateSelectizeInput(session, "cluster_filter", choices = cluster_choices, selected = character(0), server = FALSE)
      updateSelectizeInput(session, "gene_filter", choices = gene_choices, selected = character(0), server = FALSE)
    }, ignoreInit = FALSE)

    filtered_deg <- reactive({
      df <- deg_table_raw()

      filter_precomputed_deg_table(
        df = df,
        cluster_filter = input$cluster_filter %||% character(0),
        gene_filter = input$gene_filter %||% character(0),
        max_p_adj = input$max_p_adj,
        min_abs_log2fc = input$min_abs_log2fc,
        min_pct_1 = input$min_pct_1,
        min_pct_2 = input$min_pct_2,
        sort_by = input$sort_by,
        sort_desc = input$sort_desc
      )
    })

    output$deg_status <- renderText({
      path <- deg_workbook_path()
      validate(need(!is.null(path) && file.exists(path), "No precomputed DEG workbook is mapped for this Seurat object yet."))

      paste0("Using precomputed workbook: ", basename(path))
    })

    output$deg_table <- renderDataTable({
      df <- filtered_deg()
      validate(need(nrow(df) > 0, "No DEG rows remain after applying the selected filters."))
      format_deg_table_for_display(df)
    }, options = list(pageLength = 25, scrollX = TRUE))

    output$quadrant_status <- renderText({
      cache <- quadrant_cache()
      validate(need(!is.null(cache), "Quadrant scatter cache is not available for this object yet. Build it from the full Atlas object first."))
      sample_col <- default_quadrant_sample_col(loaded()$obj)
      validate(need(!is.null(sample_col), "Metadata column 'sampleID' is required for the quadrant scatter."))
      denom_tbl <- quadrant_sample_denominators(
        obj = loaded()$obj,
        sample_col = sample_col,
        group_col = cache$group_col %||% "final_group"
      )
      denom_text <- paste(
        paste0(denom_tbl$quadrant, ": n=", denom_tbl$n_samples),
        collapse = " | "
      )
      paste0(
        "Logfold changes of genes in malignant vs primary epithelial cells. ",
        "Dot size by percentage of samples expressing minimum threshold value (default 0.1). ",
        "Denominators: ",
        denom_text,
        "."
      )
    })

    output$quadrant_plot <- renderPlot({
      req(quadrant_results())
      plot_quadrant_deg(quadrant_results())
    })

    output$quadrant_table <- renderDataTable({
      req(quadrant_results())
      quadrant_results()$plot_df |>
        dplyr::arrange(dplyr::desc(distance))
    }, options = list(pageLength = 20, scrollX = TRUE))

    output$download_deg <- downloadHandler(
      filename = function() {
        paste0("precomputed_deg_filtered.csv")
      },
      content = function(file) {
        write.csv(filtered_deg(), file, row.names = FALSE)
      }
    )
  })
}
