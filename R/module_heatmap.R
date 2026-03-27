heatmap_ui <- function(id) {
  ns <- NS(id)

  fluidPage(
    fluidRow(
      column(
        width = 3,
        selectizeInput(
          ns("genes"),
          "Genes for average expression heatmap",
          choices = NULL,
          multiple = TRUE,
          options = list(
            placeholder = "Select genes",
            maxOptions = 100000
          )
        ),
        selectInput(ns("display_by"), "Display columns by", choices = character(0)),
        uiOutput(ns("filter_controls")),
        tags$div(
          class = "small-note heatmap-note",
          "This panel uses precomputed sample-by-celltype averages when a cache is available, then filters and regroups that cache inside the app."
        )
      ),
      column(
        width = 9,
        card(
          card_header("Average Expression Heatmap"),
          plotOutput(ns("heatmap"), height = 520)
        ),
        br(),
        card(
          card_header("Average Expression Table"),
          dataTableOutput(ns("avg_table"))
        )
      )
    )
  )
}

heatmap_server <- function(id, loaded) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    cache_data <- reactive({
      req(loaded()$source_path)
      if (!average_heatmap_cache_exists(loaded()$source_path)) {
        return(NULL)
      }
      read_average_heatmap_cache(loaded()$source_path)
    })

    filter_values <- reactive({
      cache <- cache_data()
      if (is.null(cache)) {
        return(list())
      }
      fields <- available_heatmap_filter_fields(cache)
      stats::setNames(
        lapply(fields, function(field) input[[paste0("filter_", field)]] %||% character(0)),
        fields
      )
    })

    avg_data <- reactive({
      req(loaded()$obj)
      req(input$display_by)
      validate(need(length(input$genes) > 0, "Choose at least one gene."))

      build_average_heatmap_data(
        genes = input$genes,
        source_path = loaded()$source_path,
        obj = loaded()$obj,
        display_by = input$display_by,
        filters = filter_values()
      )
    })

    observeEvent(loaded()$obj, {
      req(loaded()$obj)

      available_gene_choices <- if (!is.null(cache_data())) sort(rownames(cache_data()$avg_mat)) else available_features(loaded()$obj)
      updateSelectizeInput(
        session,
        "genes",
        choices = available_gene_choices,
        selected = default_heatmap_genes(available_gene_choices),
        server = FALSE
      )
      if (!is.null(cache_data())) {
        updateSelectInput(
          session,
          "display_by",
          choices = build_heatmap_display_choices(cache_data()),
          selected = default_heatmap_display_field(cache_data())
        )
      } else {
        updateSelectInput(
          session,
          "display_by",
          choices = c("final_celltype x final_group (live)" = "sample_celltype"),
          selected = "sample_celltype"
        )
      }
    })

    output$filter_controls <- renderUI({
      cache <- cache_data()
      if (is.null(cache)) {
        return(tags$div(class = "small-note", "No average-expression cache found yet for this object. Build the cache to enable fast metadata filtering."))
      }

      fields <- available_heatmap_filter_fields(cache)
      tagList(
        lapply(fields, function(field) {
          field_choices <- sort(unique(as.character(cache$meta[[field]])))
          selectizeInput(
            ns(paste0("filter_", field)),
            paste0(heatmap_filter_labels[[field]], " filter"),
            choices = field_choices,
            selected = default_heatmap_filter_values(cache, field),
            multiple = TRUE,
            options = list(placeholder = "Leave empty to keep all")
          )
        })
      )
    })

    output$heatmap <- renderPlot({
      req(avg_data())
      plot_average_heatmap_from_script(avg_data())
    })

    output$avg_table <- renderDataTable({
      req(avg_data())
      avg_data()$long
    }, options = list(pageLength = 12, scrollX = TRUE))
  })
}
