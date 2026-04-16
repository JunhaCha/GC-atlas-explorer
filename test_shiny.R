#!/usr/bin/env Rscript

script_path <- normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]),
  winslash = "/",
  mustWork = FALSE
)
repo_dir <- dirname(script_path)

app_data_dir <- function() {
  env_dir <- Sys.getenv("GC_APP_DATA_DIR", unset = "")
  if (nzchar(env_dir)) {
    return(normalizePath(env_dir, winslash = "/", mustWork = FALSE))
  }

  local_data_dir <- file.path(repo_dir, "data")
  if (dir.exists(local_data_dir)) {
    return(normalizePath(local_data_dir, winslash = "/", mustWork = FALSE))
  }

  normalizePath(dirname(repo_dir), winslash = "/", mustWork = FALSE)
}

message("Launching app from: ", repo_dir)
message("Using data dir: ", app_data_dir())

Sys.setenv(GC_APP_DATA_DIR = app_data_dir())
setwd(repo_dir)
shiny::runApp(".", launch.browser = TRUE)
