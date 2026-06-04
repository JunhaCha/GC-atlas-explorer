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

atlas_filenames <- c(
  "seurat_merged_TME_malignant_final_umap_app_slim.rds"
)

atlas_paths <- file.path(app_data_dir(), atlas_filenames)
build_script <- file.path(repo_dir, "build_quadrant_cache_for_app.R")
force_rebuild <- "--force" %in% commandArgs(trailingOnly = TRUE)
rscript_bin <- Sys.which("Rscript")
if (!nzchar(rscript_bin)) {
  stop("Rscript was not found on PATH.")
}

for (input_path in atlas_paths) {
  if (!file.exists(input_path)) {
    message("Skipping missing slim object: ", input_path)
    next
  }

  output_path <- sub("\\.rds$", "_quadrant_cache.rds", input_path, ignore.case = TRUE)
  if (file.exists(output_path) && !force_rebuild) {
    message("Quadrant cache already exists, skipping: ", output_path)
    next
  }

  message("")
  message("==== Building Quadrant Cache ====")
  message("Input:  ", input_path)
  message("Output: ", output_path)

  status <- system2(
    command = rscript_bin,
    args = c(build_script, "--input", input_path, "--output", output_path),
    env = c("OMP_NUM_THREADS=1", "R_MAX_VSIZE=100Gb"),
    stdout = "",
    stderr = ""
  )

  if (!identical(status, 0L)) {
    stop("Quadrant cache build failed for: ", input_path)
  }
}

message("")
message("Finished checking/building quadrant caches.")
