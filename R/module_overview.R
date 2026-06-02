overview_ui <- function(id) {
  ns <- NS(id)

  fluidPage(
    fluidRow(
      column(
        width = 5,
        card(
          card_header("Object Summary"),
          tableOutput(ns("summary_table"))
        )
      ),
      column(
        width = 7,
        card(
          card_header("Metadata Preview"),
          DT::DTOutput(ns("metadata_table"))
        )
      )
    )
  )
}

overview_server <- function(id, loaded) {
  moduleServer(id, function(input, output, session) {
    output$summary_table <- renderTable({
      req(loaded()$obj)
      object_summary(loaded()$obj, loaded()$path, loaded()$group_var)
    }, rownames = FALSE, colnames = FALSE)

    output$metadata_table <- DT::renderDT({
      req(loaded()$obj)
      metadata_preview(loaded()$obj)
    },
    server = TRUE,
    options = datatable_simple_pager_options(
      page_length = 12,
      extra = list(
        searching = FALSE,
        orderClasses = TRUE
      )
    ))
  })
}
