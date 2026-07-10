xenium_overview_ui <- function(id) {
  ns <- NS(id)

  fluidPage(
    fluidRow(
      column(
        width = 4,
        card(
          card_header("Sample Summary"),
          tags$p(class = "small-note", xenium_direct_neighborhood_definition()),
          tableOutput(ns("summary_table"))
        )
      ),
      column(
        width = 4,
        card(
          card_header("Cell Type Counts"),
          DT::DTOutput(ns("celltype_table"))
        )
      ),
      column(
        width = 4,
        card(
          card_header("Direct Neighborhood Counts"),
          DT::DTOutput(ns("neighborhood_table"))
        )
      )
    ),
    fluidRow(
      column(
        width = 6,
        card(
          card_header("Cells by Cell Type"),
          plotOutput(ns("celltype_plot"), height = 360)
        )
      ),
      column(
        width = 6,
        card(
          card_header("Cells by Direct Neighborhood"),
          plotOutput(ns("neighborhood_plot"), height = 360)
        )
      )
    )
  )
}

xenium_overview_server <- function(id, loaded) {
  moduleServer(id, function(input, output, session) {
    bundle <- reactive({
      req(loaded()$bundle)
      loaded()$bundle
    })

    output$summary_table <- renderTable({
      xenium_sample_summary(bundle())
    }, rownames = FALSE, colnames = FALSE)

    output$celltype_table <- DT::renderDT({
      xenium_count_table(bundle()$cells, "celltype")
    }, server = TRUE, options = datatable_simple_pager_options(page_length = 10))

    output$neighborhood_table <- DT::renderDT({
      xenium_count_table(bundle()$cells, "neighborhood")
    }, server = TRUE, options = datatable_simple_pager_options(page_length = 10))

    output$celltype_plot <- renderPlot({
      xenium_overview_barplot(bundle()$cells, "celltype", "Cells by cell type")
    })

    output$neighborhood_plot <- renderPlot({
      xenium_overview_barplot(bundle()$cells, "neighborhood", "Cells by direct neighborhood")
    })
  })
}
