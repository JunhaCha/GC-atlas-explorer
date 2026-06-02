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
          dataTableOutput(ns("metadata_table"))
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

    output$metadata_table <- renderDataTable({
      req(loaded()$obj)
      metadata_preview(loaded()$obj)
    }, options = list(
      scrollX = TRUE,
      paging = FALSE,
      searching = FALSE,
      info = FALSE,
      dom = "t",
      orderClasses = TRUE
    ))
  })
}
