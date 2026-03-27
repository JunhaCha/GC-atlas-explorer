source("R/helpers.R")
source("R/module_data_loader.R")
source("R/module_overview.R")
source("R/module_umap.R")
source("R/module_gene_explorer.R")
source("R/module_heatmap.R")
source("R/module_deg.R")

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(Seurat)
})

atlas_object_choices <- atlas_object_choices_default()

xenium_object_choices <- xenium_object_choices_default()

dataset_workspace_ui <- function(prefix, title, default_path = "", note = NULL, object_choices = NULL) {
  layout_sidebar(
    sidebar = sidebar(
      width = 340,
      data_loader_ui(
        id = paste0(prefix, "_loader"),
        title = paste(title, "Data"),
        default_path = default_path,
        note = note,
        object_choices = object_choices
      )
    ),
    navset_tab(
      id = paste0(prefix, "_tabs"),
      nav_panel("Overview", overview_ui(paste0(prefix, "_overview"))),
      nav_panel("UMAP", umap_ui(paste0(prefix, "_umap"))),
      nav_panel("Gene Explorer", gene_explorer_ui(paste0(prefix, "_gene_explorer"))),
      nav_panel("Sample Explorer", heatmap_ui(paste0(prefix, "_heatmap"))),
      nav_panel("DEG", deg_ui(paste0(prefix, "_deg")))
    )
  )
}

ui <- page_navbar(
  title = "Gastric Cancer single-cell Explorer",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  header = tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
  ),
  nav_panel(
    "scRNA-seq Atlas",
    dataset_workspace_ui(
      prefix = "atlas",
      title = "scRNA-seq Atlas",
      object_choices = atlas_object_choices
    )
  ),
  nav_panel(
    "Xenium Spatial",
    dataset_workspace_ui(
      prefix = "xenium",
      title = "Xenium Spatial",
      note = "Xenium objects will be wired into this selector once they are copied into Playground.",
      object_choices = xenium_object_choices
    )
  )
)

server <- function(input, output, session) {
  atlas_loaded <- data_loader_server("atlas_loader", object_choices = atlas_object_choices)
  xenium_loaded <- data_loader_server("xenium_loader", object_choices = xenium_object_choices)

  overview_server("atlas_overview", atlas_loaded)
  umap_server("atlas_umap", atlas_loaded)
  gene_explorer_server("atlas_gene_explorer", atlas_loaded)
  heatmap_server("atlas_heatmap", atlas_loaded)
  deg_server("atlas_deg", atlas_loaded)

  overview_server("xenium_overview", xenium_loaded)
  umap_server("xenium_umap", xenium_loaded)
  gene_explorer_server("xenium_gene_explorer", xenium_loaded)
  heatmap_server("xenium_heatmap", xenium_loaded)
  deg_server("xenium_deg", xenium_loaded)
}

shinyApp(ui, server)
