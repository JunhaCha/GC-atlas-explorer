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
            placeholder = "Search genes",
            maxOptions = 100000,
            maxItems = 4
          )
        ),
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
          plotOutput(
            ns("feature_plot"),
            height = 320,
            brush = brushOpts(ns("feature_brush"), direction = "xy", resetOnNew = TRUE)
          ),
          div(
            class = "mt-2",
            actionButton(ns("reset_feature_zoom"), "Reset zoom", class = "btn-outline-secondary btn-sm")
          )
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

violin_choice_blacklist <- c(
  "age",
  "sampleID",
  "patientID",
  "nCount_RNA",
  "nFeature_RNA",
  "percent.mito",
  "nCount_SCT",
  "nFeature_SCT"
)

filtered_violin_choices <- function(obj) {
  cols <- colnames(obj[[]])
  cols[!cols %in% violin_choice_blacklist]
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
    validate(need("final_group" %in% colnames(out[[]]), "final_group is required for the default violin grouping."))
    out$.violin_pathology_group <- derive_violin_pathology_group(out$final_group)
  }

  if (identical(group_var, ".violin_gc_status") || identical(split_var, ".violin_gc_status")) {
    validate(need("final_group" %in% colnames(out[[]]), "final_group is required for the default violin split."))
    out$.violin_gc_status <- derive_violin_gc_status(out$final_group)
  }

  out
}

gene_explorer_server <- function(id, loaded) {
  moduleServer(id, function(input, output, session) {
    feature_zoom <- reactiveValues(x = NULL, y = NULL)

    observeEvent(loaded()$obj, {
      req(loaded()$obj)
      feature_zoom$x <- NULL
      feature_zoom$y <- NULL
      violin_choices <- filtered_violin_choices(loaded()$obj)

      updateSelectizeInput(
        session,
        "genes",
        choices = available_features(loaded()$obj),
        selected = {
          features <- available_features(loaded()$obj)
          if ("A1BG" %in% features) "A1BG" else head(features, 1)
        },
        server = FALSE
      )
      updateSelectInput(
        session,
        "group_by",
        choices = c(
          "Pathological Group (from final_group)" = ".violin_pathology_group",
          stats::setNames(violin_choices, violin_choices)
        ),
        selected = {
          if ("final_group" %in% colnames(loaded()$obj[[]])) {
            ".violin_pathology_group"
          } else {
            preferred <- c("rev_pathological_subtype")
            hits <- preferred[preferred %in% violin_choices]
            if (length(hits) > 0) hits[[1]] else loaded()$group_var
          }
        }
      )
      updateSelectInput(
        session,
        "split_by",
        choices = c(
          "None" = "",
          "Normal vs matched GC (from final_group)" = ".violin_gc_status",
          stats::setNames(violin_choices, violin_choices)
        ),
        selected = if ("final_group" %in% colnames(loaded()$obj[[]])) ".violin_gc_status" else if ("rev_condition" %in% violin_choices) "rev_condition" else ""
      )
      if ("final_celltype" %in% colnames(loaded()$obj[[]])) {
        updateSelectizeInput(
          session,
          "violin_celltypes",
          choices = sort(unique(as.character(loaded()$obj$final_celltype))),
          selected = character(0),
          server = FALSE
        )
      } else {
        updateSelectizeInput(
          session,
          "violin_celltypes",
          choices = character(0),
          selected = character(0),
          server = FALSE
        )
      }
    })

    observeEvent(input$genes, {
      feature_zoom$x <- NULL
      feature_zoom$y <- NULL
    }, ignoreInit = TRUE)

    observeEvent(input$feature_brush, {
      brush <- input$feature_brush
      req(brush)
      if (is.null(brush$xmin) || is.null(brush$xmax) || is.null(brush$ymin) || is.null(brush$ymax)) {
        return(invisible(NULL))
      }
      feature_zoom$x <- c(brush$xmin, brush$xmax)
      feature_zoom$y <- c(brush$ymin, brush$ymax)
    })

    observeEvent(input$reset_feature_zoom, {
      feature_zoom$x <- NULL
      feature_zoom$y <- NULL
    }, ignoreInit = TRUE)

    output$feature_plot <- renderPlot({
      req(loaded()$obj, input$genes)
      validate(need(length(input$genes) > 0, "Choose at least one gene."))
      validate(need(length(input$genes) <= 4, "Choose up to 4 genes."))

      plot_feature_expression(
        obj = loaded()$obj,
        features = input$genes[seq_len(min(4, length(input$genes)))],
        xlim = feature_zoom$x,
        ylim = feature_zoom$y
      )
    })

    output$violin_plot <- renderPlot({
      req(loaded()$obj, input$genes, input$group_by)
      validate(need(length(input$genes) <= 4, "Choose up to 4 genes."))
      split_var <- if (nzchar(input$split_by %||% "")) input$split_by else NULL
      obj_to_plot <- loaded()$obj
      selected_celltypes <- input$violin_celltypes %||% character(0)
      if (length(selected_celltypes) > 0) {
        validate(need("final_celltype" %in% colnames(obj_to_plot[[]]), "final_celltype column is not available for violin filtering."))
        keep_cells <- rownames(obj_to_plot[[]])[as.character(obj_to_plot$final_celltype) %in% selected_celltypes]
        validate(need(length(keep_cells) > 0, "No cells remain after applying the selected cell-type filter."))
        obj_to_plot <- subset(obj_to_plot, cells = keep_cells)
      }
      obj_to_plot <- materialize_violin_grouping(obj_to_plot, input$group_by, split_var)
      validate(need(input$group_by %in% colnames(obj_to_plot[[]]), "Selected grouping column is not available."))
      if (!is.null(split_var)) {
        validate(need(split_var %in% colnames(obj_to_plot[[]]), "Selected split column is not available."))
      }

      plot_violin_expression(
        obj = obj_to_plot,
        features = input$genes[seq_len(min(4, length(input$genes)))],
        group_var = input$group_by,
        split_var = split_var
      )
    })
  })
}
