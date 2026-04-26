xenium_marker_ui <- function(id) {
  ns <- NS(id)

  fluidPage(
    fluidRow(
      column(
        width = 3,
        card(
          card_header("Marker Overlay Controls"),
          selectizeInput(
            ns("gene"),
            "Gene",
            choices = character(0),
            selected = character(0),
            options = list(placeholder = "Choose a Xenium gene")
          ),
          sliderInput(ns("point_size"), "Point size", min = 0.1, max = 2.5, value = 0.5, step = 0.1),
          sliderInput(ns("alpha"), "Point alpha", min = 0.1, max = 1, value = 0.8, step = 0.05),
          tags$div(
            class = "small-note",
            "Loads one Xenium gene at a time from the app-ready sparse export."
          )
        )
      ),
      column(
        width = 9,
        card(
          card_header("Marker Overlay"),
          plotOutput(ns("marker_plot"), height = 780)
        )
      )
    )
  )
}

xenium_marker_server <- function(id, loaded) {
  moduleServer(id, function(input, output, session) {
    chosen_gene <- reactiveVal(NULL)

    bundle <- reactive({
      req(loaded()$bundle)
      loaded()$bundle
    })

    available_genes <- reactive({
      req(bundle()$sample_id)
      tryCatch(
        xenium_marker_available_genes(bundle()$sample_id),
        error = function(e) character(0)
      )
    })

    observeEvent(input$gene, {
      if (is.null(input$gene)) {
        return()
      }
      chosen_gene(input$gene)
    }, ignoreInit = FALSE)

    observeEvent(bundle(), {
      genes <- available_genes()
      remembered_gene <- chosen_gene() %||% character(0)
      selected_gene <- if (
        length(remembered_gene) == 1 &&
        nzchar(remembered_gene) &&
        remembered_gene %in% genes
      ) {
        remembered_gene
      } else {
        character(0)
      }

      updateSelectizeInput(
        session,
        "gene",
        choices = genes,
        selected = selected_gene,
        server = FALSE
      )
    }, ignoreInit = FALSE)

    output$marker_plot <- renderPlot({
      validate(need(length(available_genes()) > 0, paste0(
        "Marker assets are not available for sample ",
        bundle()$sample_id,
        ". Rebuild them with ",
        app_script_path("build_xenium_marker_assets.R"),
        "."
      )))
      validate(need(
        !is.null(input$gene) && nzchar(input$gene),
        "Choose a gene to display the marker overlay."
      ))
      plot_xenium_marker_overlay(
        bundle = bundle(),
        gene = input$gene,
        point_size = input$point_size %||% 0.5,
        alpha = input$alpha %||% 0.8
      )
    }, res = 110)
  })
}
