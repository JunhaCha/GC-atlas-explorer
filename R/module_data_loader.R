data_loader_ui <- function(
  id,
  title = "Data",
  default_path = "",
  note = NULL,
  object_choices = NULL
) {
  ns <- NS(id)
  use_selector <- !is.null(object_choices)

  tagList(
    h4(title),
    if (use_selector) {
      selectInput(
        ns("object_choice"),
        "Choose Seurat object",
        choices = object_choices,
        selected = if (length(object_choices) > 0) object_choices[[1]] else character(0)
      )
    } else {
      textInput(
        ns("rds_path"),
        "Local path to Seurat .rds",
        value = default_path,
        placeholder = "/absolute/path/to/object.rds"
      )
    },
    actionButton(
      ns("load"),
      if (use_selector) "Use selected object" else "Load object",
      class = "btn-primary"
    ),
    tags$hr(),
    selectInput(ns("group_var"), "Primary group / identity column", choices = character(0)),
    selectInput(ns("assay"), "Default assay", choices = character(0)),
    tags$hr(),
    strong("Status"),
    textOutput(ns("status")),
    if (!is.null(note) && nzchar(note)) {
      tags$div(class = "small-note", note)
    }
  )
}

data_loader_server <- function(id, object_choices = NULL) {
  moduleServer(id, function(input, output, session) {
    values <- reactiveValues(
      obj = NULL,
      path = NULL,
      source_path = NULL,
      status = "Load a Seurat .rds file to begin."
    )

    load_selected_object <- function() {
      path_to_use <- NULL
      choice_map <- object_choices %||% NULL

      if (!is.null(choice_map)) {
        if (!is.null(input$object_choice) && nzchar(input$object_choice)) {
          path_to_use <- input$object_choice
          values$path <- names(choice_map)[match(input$object_choice, unname(choice_map))] %||% input$object_choice
        }
      } else if (nzchar(trimws(input$rds_path %||% ""))) {
        path_to_use <- trimws(input$rds_path)
        values$path <- trimws(input$rds_path)
      }

      if (is.null(path_to_use)) {
        values$status <- "Choose a Seurat object."
        return(invisible(NULL))
      }

      tryCatch({
        obj <- safe_seurat_read(path_to_use)

        updateSelectInput(session, "group_var", choices = colnames(obj[[]]), selected = default_group_var(obj))
        updateSelectInput(session, "assay", choices = available_assays(obj), selected = DefaultAssay(obj))

        values$obj <- obj
        values$source_path <- path_to_use
        values$status <- paste0(
          "Loaded object with ",
          format(ncol(obj), big.mark = ","),
          " cells and ",
          format(nrow(obj), big.mark = ","),
          " features."
        )
      }, error = function(e) {
        values$obj <- NULL
        values$source_path <- path_to_use
        values$status <- conditionMessage(e)
      })
    }

    observeEvent(input$load, {
      load_selected_object()
    }, ignoreInit = TRUE)

    observeEvent(input$object_choice, {
      if (!is.null(object_choices) && !is.null(input$object_choice) && nzchar(input$object_choice)) {
        load_selected_object()
      }
    }, ignoreInit = TRUE)

    observeEvent(list(input$group_var, input$assay), {
      req(values$obj)

      if (!is.null(input$assay) && nzchar(input$assay)) {
        DefaultAssay(values$obj) <- input$assay
      }

      if (!is.null(input$group_var) && nzchar(input$group_var)) {
        Idents(values$obj) <- values$obj[[input$group_var, drop = TRUE]]
      }
    })

    output$status <- renderText(values$status)

    reactive({
      list(
        obj = values$obj,
        path = values$path,
        source_path = values$source_path,
        group_var = input$group_var,
        assay = input$assay
      )
    })
  })
}
