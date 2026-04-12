xenium_loader_ui <- function(id, title = "Xenium Spatial Data") {
  ns <- NS(id)

  choices <- tryCatch(
    xenium_sample_choices(),
    error = function(e) c("No Xenium samples found" = "")
  )

  tagList(
    h4(title),
    selectInput(
      ns("sample_id"),
      "Choose sample",
      choices = choices,
      selected = if (length(choices) > 0) choices[[1]] else character(0)
    ),
    tags$hr(),
    strong("Status"),
    textOutput(ns("status")),
    tags$div(
      class = "small-note",
      "Loads only the selected sample image, cell centroids, and precomputed neighborhood summaries."
    )
  )
}

xenium_loader_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    values <- reactiveValues(
      bundle = NULL,
      status = "Choose a Xenium sample to begin."
    )

    manifest <- reactive({
      xenium_sample_manifest()
    })

    observeEvent(input$sample_id, {
      req(input$sample_id, nzchar(input$sample_id))

      tryCatch({
        bundle <- xenium_sample_bundle(input$sample_id, manifest())
        values$bundle <- bundle
        values$status <- paste0(
          "Loaded sample ",
          input$sample_id,
          " with ",
          format(nrow(bundle$cells), big.mark = ","),
          " cells."
        )
      }, error = function(e) {
        values$bundle <- NULL
        values$status <- conditionMessage(e)
      })
    }, ignoreInit = FALSE)

    output$status <- renderText(values$status)

    reactive({
      list(
        sample_id = input$sample_id,
        manifest = tryCatch(manifest(), error = function(e) NULL),
        bundle = values$bundle,
        status = values$status
      )
    })
  })
}
