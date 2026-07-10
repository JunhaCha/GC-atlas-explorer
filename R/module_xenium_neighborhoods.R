xenium_neighborhoods_ui <- function(id) {
  ns <- NS(id)

  fluidPage(
    tags$p(class = "small-note", xenium_direct_neighborhood_definition()),
    fluidRow(
      column(
        width = 5,
        card(
          card_header("Direct Neighborhood Abundance"),
          plotOutput(ns("abundance_plot"), height = 420)
        )
      ),
      column(
        width = 7,
        card(
          card_header("Direct Neighborhood Composition"),
          plotOutput(ns("composition_plot"), height = 420)
        )
      )
    ),
    fluidRow(
      column(
        width = 5,
        card(
          card_header("Direct Neighborhood Summary Table"),
          dataTableOutput(ns("summary_table"))
        )
      ),
      column(
        width = 7,
        card(
          card_header("Direct Neighborhood Composition Table"),
          dataTableOutput(ns("composition_table"))
        )
      )
    )
  )
}

xenium_neighborhoods_server <- function(id, loaded) {
  moduleServer(id, function(input, output, session) {
    bundle <- reactive({
      req(loaded()$bundle)
      loaded()$bundle
    })

    output$abundance_plot <- renderPlot({
      plot_xenium_neighborhood_abundance(bundle()$neighborhood_summary)
    })

    output$composition_plot <- renderPlot({
      plot_xenium_neighborhood_composition(
        bundle()$neighborhood_composition,
        bundle()$neighborhood_summary
      )
    })

    output$summary_table <- renderDataTable({
      bundle()$neighborhood_summary |>
        dplyr::mutate(fraction = sprintf("%.1f%%", fraction * 100))
    }, options = datatable_simple_pager_options(page_length = 10))

    output$composition_table <- renderDataTable({
      bundle()$neighborhood_composition |>
        dplyr::mutate(fraction = sprintf("%.1f%%", fraction * 100))
    }, options = datatable_simple_pager_options(page_length = 12))
  })
}
