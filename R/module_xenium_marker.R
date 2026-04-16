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

    observeEvent(bundle(), {
      genes <- available_genes()

      default_gene <- if ("MIF" %in% genes) {
        "MIF"
      } else if (length(genes) > 0) {
        genes[[1]]
      } else {
        character(0)
      }

      updateSelectizeInput(
        session,
        "gene",
        choices = genes,
        selected = default_gene,
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
      req(input$gene, nzchar(input$gene))
      plot_xenium_marker_overlay(
        bundle = bundle(),
        gene = input$gene,
        point_size = input$point_size %||% 0.5,
        alpha = input$alpha %||% 0.8
      )
    }, res = 110)
  })
}
