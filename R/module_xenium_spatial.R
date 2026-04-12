xenium_spatial_ui <- function(id) {
  ns <- NS(id)

  fluidPage(
    fluidRow(
      column(
        width = 3,
        card(
          card_header("Spatial Map Controls"),
          selectInput(ns("color_by"), "Color cells by", choices = character(0)),
          checkboxInput(ns("highlight_mode"), "Highlight selected cell types", value = FALSE),
          selectInput(ns("highlight_a"), "Highlight cell type 1", choices = character(0), selected = ""),
          selectInput(ns("highlight_b"), "Highlight cell type 2", choices = character(0), selected = ""),
          sliderInput(ns("point_size"), "Point size", min = 0.1, max = 2.5, value = 0.5, step = 0.1),
          sliderInput(ns("alpha"), "Point alpha", min = 0.1, max = 1, value = 0.8, step = 0.05)
        )
      ),
      column(
        width = 9,
        card(
          card_header("Spatial Map"),
          plotOutput(ns("spatial_plot"), height = 780)
        )
      )
    )
  )
}

xenium_spatial_server <- function(id, loaded) {
  moduleServer(id, function(input, output, session) {
    bundle <- reactive({
      req(loaded()$bundle)
      loaded()$bundle
    })

    observeEvent(bundle(), {
      available_fields <- c("celltype", "subtype", "tumor_normal", "neighborhood")
      available_fields <- available_fields[available_fields %in% colnames(bundle()$cells)]
      choices <- stats::setNames(available_fields, xenium_color_labels[available_fields] %||% available_fields)
      celltype_choices <- sort(unique(as.character(bundle()$cells$celltype)))
      celltype_choices <- celltype_choices[!is.na(celltype_choices) & nzchar(celltype_choices)]
      highlight_choices <- c(stats::setNames("", ""), stats::setNames(celltype_choices, celltype_choices))

      updateSelectInput(
        session,
        "color_by",
        choices = choices,
        selected = if ("celltype" %in% available_fields) "celltype" else available_fields[[1]]
      )
      updateSelectInput(session, "highlight_a", choices = highlight_choices, selected = "")
      updateSelectInput(session, "highlight_b", choices = highlight_choices, selected = "")
    }, ignoreInit = FALSE)

    output$spatial_plot <- renderPlot({
      highlight_values <- character(0)
      if (isTRUE(input$highlight_mode)) {
        highlight_values <- c(input$highlight_a %||% "", input$highlight_b %||% "")
      }

      plot_xenium_spatial_map(
        bundle = bundle(),
        color_by = input$color_by %||% "celltype",
        point_size = input$point_size %||% 0.5,
        alpha = input$alpha %||% 0.8,
        highlight_celltypes = highlight_values
      )
    }, res = 110)
  })
}
