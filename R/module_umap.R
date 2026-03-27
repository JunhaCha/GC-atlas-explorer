umap_ui <- function(id) {
  ns <- NS(id)

  fluidPage(
    fluidRow(
      column(
        width = 3,
        tags$p(
          class = "small-note",
          "UMAP is drawn from the HarmonyUMAP reduction stored in the Seurat object."
        ),
        selectInput(
          ns("color_mode"),
          "Color by",
          choices = c(
            "Cell Type" = "final_celltype",
            "Pathological Group" = "rev_pathological_subtype",
            "Primary/Normal" = "rev_condition",
            "MSS/MSI" = "rev_molecular_subtype",
            "Cohort" = "dataset",
            "Sample ID" = "sampleID",
            "Patient ID" = "patientID",
            "TNM" = "TNM",
            "Phase" = "Phase"
          ),
          selected = "final_celltype"
        ),
        checkboxInput(ns("label"), "Show labels", TRUE),
        sliderInput(ns("pt_size"), "Point size", min = 0.1, max = 2.5, value = 1.5, step = 0.1)
      ),
      column(
        width = 9,
        card(
          card_header("Embedding Plot"),
          plotOutput(ns("umap_plot"), height = 620)
        )
      )
    )
  )
}

umap_server <- function(id, loaded) {
  moduleServer(id, function(input, output, session) {
    output$umap_plot <- renderPlot({
      req(loaded()$obj, input$color_mode)
      plot_harmony_umap(
        obj = loaded()$obj,
        color_mode = input$color_mode,
        label = input$label,
        pt_size = input$pt_size
      )
    })
  })
}
