gene_explorer_ui <- function(id) {
  ns <- NS(id)

  fluidPage(
    fluidRow(
      column(
        width = 3,
        selectizeInput(
          ns("genes"),
          "Genes",
          choices = NULL,
          multiple = TRUE,
          options = list(
            placeholder = "Select up to 4 genes",
            maxOptions = 200,
            loadThrottle = 150,
            minLength = 2,
            maxItems = 4
          )
        ),
        actionButton(
          ns("plot_genes"),
          "Plot selected gene(s)",
          class = "btn-primary w-100"
        ),
        br(),
        selectInput(ns("group_by"), "Violin x-axis grouping", choices = character(0)),
        selectInput(ns("split_by"), "Split violin within group", choices = character(0)),
        selectizeInput(
          ns("violin_celltypes"),
          "Filter violin to cell types",
          choices = NULL,
          multiple = TRUE,
          options = list(placeholder = "Leave empty to use all cell types")
        )
      ),
      column(
        width = 9,
        card(
          card_header("Feature Plot"),
          plotOutput(ns("feature_plot"), height = 320)
        ),
        br(),
        card(
          card_header("Violin Plot"),
          plotOutput(ns("violin_plot"), height = 320)
        )
      )
    )
  )
}

violin_choice_allow_patterns <- c(
  "^final_group$",
  "^final_celltype$",
  "^sample_celltype$",
  "^rev_condition$",
  "^rev_pathological_subtype$",
  "^rev_molecular_subtype$",
  "^rev_stage$",
  "^dataset$",
  "^cohort$",
  "^patientID$",
  "^sampleID$",
  "^TNM$",
  "^Pathologic_[TNM]$",
  "^sex$",
  "^gender$",
  "^tumor_normal$",
  "cell.?type",
  "patholog",
  "molecular",
  "condition",
  "stage"
)

violin_choice_block_patterns <- c(
  "^X(\\.\\d+)?$",
  "^age$",
  "^orig\\.ident$",
  "^seurat_clusters$",
  "snn",
  "^nCount_",
  "^nFeature_",
  "^percent",
  "mito",
  "ribosomal",
  "doublet",
  "barcode",
  "cell[_\\. ]?id"
)

filtered_violin_choices <- function(obj) {
  md <- atlas_meta(obj)
  cols <- colnames(md)
  allowed <- vapply(
    cols,
    function(col) any(grepl(paste(violin_choice_allow_patterns, collapse = "|"), col, ignore.case = TRUE)),
    logical(1)
  )
  blocked <- vapply(
    cols,
    function(col) any(grepl(paste(violin_choice_block_patterns, collapse = "|"), col, ignore.case = TRUE)),
    logical(1)
  )
  numeric_cols <- vapply(md, is.numeric, logical(1))

  cols[allowed & !blocked & !numeric_cols]
}

label_violin_choices <- function(cols) {
  stats::setNames(
    cols,
    vapply(cols, violin_display_label, character(1), USE.NAMES = FALSE)
  )
}

derive_violin_pathology_group <- function(final_group) {
  x <- as.character(final_group)
  out <- rep(NA_character_, length(x))
  out[grepl("^Diffuse", x)] <- "Diffuse"
  out[grepl("^Intestinal", x)] <- "Intestinal"
  out[grepl("^Mixed", x)] <- "Mixed"
  out
}

derive_violin_gc_status <- function(final_group) {
  x <- as.character(final_group)
  out <- rep(NA_character_, length(x))
  out[grepl("^Diffuse.*_Normal$", x)] <- "Normal"
  out[grepl("^Intestinal.*_Normal$", x)] <- "Normal"
  out[grepl("^Mixed.*_Normal$", x)] <- "Normal"
  out[grepl("^Diffuse", x) & !grepl("_Normal$", x)] <- "Diffuse GC"
  out[grepl("^Intestinal", x) & !grepl("_Normal$", x)] <- "Intestinal GC"
  out[grepl("^Mixed", x) & !grepl("_Normal$", x)] <- "Mixed GC"
  out
}

materialize_violin_grouping <- function(obj, group_var, split_var) {
  out <- obj

  if (identical(group_var, ".violin_pathology_group") || identical(split_var, ".violin_pathology_group")) {
    validate(need("final_group" %in% colnames(atlas_meta(out)), "final_group is required for the default violin grouping."))
    out <- atlas_set_meta_column(out, ".violin_pathology_group", derive_violin_pathology_group(out$final_group))
  }

  if (identical(group_var, ".violin_gc_status") || identical(split_var, ".violin_gc_status")) {
    validate(need("final_group" %in% colnames(atlas_meta(out)), "final_group is required for the default violin split."))
    out <- atlas_set_meta_column(out, ".violin_gc_status", derive_violin_gc_status(out$final_group))
  }

  out
}

gene_explorer_server <- function(id, loaded) {
  moduleServer(id, function(input, output, session) {
    selected_genes <- reactive({
      genes <- unique(as.character(input$genes %||% character(0)))
      genes[seq_len(min(4, length(genes)))]
    }) |> shiny::debounce(250)

    plotted_genes <- reactiveVal(NULL)

    observeEvent(loaded()$obj, {
      req(loaded()$obj)
      plotted_genes(NULL)
      violin_choices <- filtered_violin_choices(loaded()$obj)
      feature_choices <- available_features(loaded()$obj, loaded()$source_path)
      has_final_group <- "final_group" %in% colnames(atlas_meta(loaded()$obj))
      synthetic_group_choices <- if (has_final_group) {
        c("Pathological Group (from final_group)" = ".violin_pathology_group")
      } else {
        character(0)
      }
      synthetic_split_choices <- if (has_final_group) {
        c("Normal vs matched GC (from final_group)" = ".violin_gc_status")
      } else {
        character(0)
      }

      updateSelectizeInput(
        session,
        "genes",
        choices = feature_choices,
        selected = character(0),
        server = TRUE
      )
      updateSelectInput(
        session,
        "group_by",
        choices = c(
          synthetic_group_choices,
          label_violin_choices(violin_choices)
        ),
        selected = {
          if (has_final_group) {
            ".violin_pathology_group"
          } else {
            preferred <- c("rev_pathological_subtype")
            hits <- preferred[preferred %in% violin_choices]
            if (length(hits) > 0) {
              hits[[1]]
            } else if (loaded()$group_var %in% violin_choices) {
              loaded()$group_var
            } else if (length(violin_choices) > 0) {
              violin_choices[[1]]
            } else {
              ""
            }
          }
        }
      )
      updateSelectInput(
        session,
        "split_by",
        choices = c(
          "None" = "",
          synthetic_split_choices,
          label_violin_choices(violin_choices)
        ),
        selected = if (has_final_group) ".violin_gc_status" else if ("rev_condition" %in% violin_choices) "rev_condition" else ""
      )
      if ("final_celltype" %in% colnames(atlas_meta(loaded()$obj))) {
        updateSelectizeInput(
          session,
          "violin_celltypes",
          choices = sort(unique(as.character(loaded()$obj$final_celltype))),
          selected = character(0),
          server = TRUE
        )
      } else {
        updateSelectizeInput(
          session,
          "violin_celltypes",
          choices = character(0),
          selected = character(0),
          server = TRUE
        )
      }
    })

    observeEvent(input$plot_genes, {
      plotted_genes(selected_genes())
    })

    output$feature_plot <- renderPlot({
      req(loaded()$obj)
      genes_to_plot <- plotted_genes()
      validate(need(!is.null(genes_to_plot), "Choose up to 4 genes and click 'Plot selected gene(s)'."))
      validate(need(length(genes_to_plot) > 0, "Choose at least one gene."))
      validate(need(length(genes_to_plot) <= 4, "Choose up to 4 genes."))

      plot_feature_expression(
        obj = loaded()$obj,
        features = genes_to_plot
      )
    })

    output$violin_plot <- renderPlot({
      req(loaded()$obj, input$group_by)
      genes_to_plot <- plotted_genes()
      validate(need(!is.null(genes_to_plot), "Choose up to 4 genes and click 'Plot selected gene(s)'."))
      validate(need(length(genes_to_plot) <= 4, "Choose up to 4 genes."))
      split_var <- if (nzchar(input$split_by %||% "")) input$split_by else NULL
      obj_to_plot <- loaded()$obj
      selected_celltypes <- input$violin_celltypes %||% character(0)
      if (length(selected_celltypes) > 0) {
        validate(need("final_celltype" %in% colnames(atlas_meta(obj_to_plot)), "final_celltype column is not available for violin filtering."))
        keep_cells <- rownames(atlas_meta(obj_to_plot))[as.character(obj_to_plot$final_celltype) %in% selected_celltypes]
        validate(need(length(keep_cells) > 0, "No cells remain after applying the selected cell-type filter."))
        obj_to_plot <- atlas_subset_cells(obj_to_plot, cells = keep_cells)
      }
      obj_to_plot <- materialize_violin_grouping(obj_to_plot, input$group_by, split_var)
      validate(need(input$group_by %in% colnames(atlas_meta(obj_to_plot)), "Selected grouping column is not available."))
      if (!is.null(split_var)) {
        validate(need(split_var %in% colnames(atlas_meta(obj_to_plot)), "Selected split column is not available."))
      }

      plot_violin_expression(
        obj = obj_to_plot,
        features = genes_to_plot,
        group_var = input$group_by,
        split_var = split_var
      )
    })
  })
}
