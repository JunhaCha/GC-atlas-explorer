source("R/helpers.R")
source("R/module_data_loader.R")
source("R/module_overview.R")
source("R/module_umap.R")
source("R/module_gene_explorer.R")
source("R/module_heatmap.R")
source("R/module_deg.R")
source("R/module_xenium_loader.R")
source("R/module_xenium_overview.R")
source("R/module_xenium_spatial.R")
source("R/module_xenium_neighborhoods.R")
source("R/module_xenium_marker.R")

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(Seurat)
})

atlas_object_choices <- atlas_object_choices_default()
preload_mode <- tolower(Sys.getenv("GC_APP_PRELOAD_ATLAS", unset = "entire"))
if (identical(preload_mode, "all")) {
  invisible(lapply(unname(atlas_object_choices), function(path) {
    try(preload_seurat_object(path), silent = TRUE)
  }))
} else if (identical(preload_mode, "none")) {
  invisible(NULL)
} else {
  try(preload_seurat_object(atlas_object_choices[["Entire Atlas"]]), silent = TRUE)
}

atlas_workspace_ui <- function(prefix, title, default_path = "", note = NULL, object_choices = NULL) {
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

xenium_workspace_ui <- function(prefix, title) {
  layout_sidebar(
    sidebar = sidebar(
      width = 340,
      xenium_loader_ui(
        id = paste0(prefix, "_loader"),
        title = paste(title, "Data")
      )
    ),
    navset_tab(
      id = paste0(prefix, "_tabs"),
      nav_panel("Sample Overview", xenium_overview_ui(paste0(prefix, "_overview"))),
      nav_panel("Neighborhood Summaries", xenium_neighborhoods_ui(paste0(prefix, "_neighborhoods"))),
      nav_panel("Spatial Map", xenium_spatial_ui(paste0(prefix, "_spatial"))),
      nav_panel("Marker Overlay", xenium_marker_ui(paste0(prefix, "_marker")))
    )
  )
}

app_brand_title <- div(
  class = "app-brand",
  span(class = "app-brand__title", "Gastric Cancer single-cell Explorer")
)

ui <- page_navbar(
  title = app_brand_title,
  window_title = "Gastric Cancer single-cell Explorer",
  theme = bs_theme(
    version = 5,
    primary = "#4E879E",
    secondary = "#9AA387",
    success = "#9AA387",
    info = "#6E9CAF",
    bg = "#F6FBFC",
    fg = "#24333B"
  ),
  header = tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
  ),
  nav_panel(
    "scRNA-seq Atlas",
    atlas_workspace_ui(
      prefix = "atlas",
      title = "scRNA-seq Atlas",
      object_choices = atlas_object_choices
    )
  ),
  nav_panel(
    "Xenium Spatial",
    xenium_workspace_ui(
      prefix = "xenium",
      title = "Xenium Spatial"
    )
  ),
  nav_spacer(),
  nav_item(
    div(
      class = "navbar-partner-brand",
      tags$img(
        src = "stanford-jrg-logo.jpeg",
        alt = "Stanford Medicine Ji Research Group",
        class = "navbar-partner-brand__logo"
      )
    )
  )
)

server <- function(input, output, session) {
  atlas_loaded <- data_loader_server("atlas_loader", object_choices = atlas_object_choices)
  xenium_loaded <- xenium_loader_server("xenium_loader")

  overview_server("atlas_overview", atlas_loaded)
  umap_server("atlas_umap", atlas_loaded)
  gene_explorer_server("atlas_gene_explorer", atlas_loaded)
  heatmap_server("atlas_heatmap", atlas_loaded)
  deg_server("atlas_deg", atlas_loaded)

  xenium_overview_server("xenium_overview", xenium_loaded)
  xenium_spatial_server("xenium_spatial", xenium_loaded)
  xenium_neighborhoods_server("xenium_neighborhoods", xenium_loaded)
  xenium_marker_server("xenium_marker", xenium_loaded)
}

shinyApp(ui, server)
