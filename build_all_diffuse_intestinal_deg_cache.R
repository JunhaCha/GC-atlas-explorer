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
  "seurat_epithelial_normal_final_final_app_slim.rds",
  "seurat_Stromal_final_app_slim.rds",
  "seurat_CD8T_final2_app_slim.rds",
  "seurat_CD4T_final2_app_slim.rds",
  "seurat_B_final2_app_slim.rds",
  "seurat_Mye_final2_app_slim.rds"
)

atlas_paths <- file.path(app_data_dir(), atlas_filenames)
build_script <- file.path(repo_dir, "build_diffuse_intestinal_deg_cache.R")
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

  output_path <- sub("\\.rds$", "_diffuse_intestinal_deg_cache.rds", input_path, ignore.case = TRUE)
  if (file.exists(output_path) && !force_rebuild) {
    message("Diffuse-vs-Intestinal DEG cache already exists, skipping: ", output_path)
    next
  }

  message("")
  message("==== Building Diffuse vs Intestinal DEG Cache ====")
  message("Input:  ", input_path)
  message("Output: ", output_path)

  status <- system2(
    command = rscript_bin,
    args = c(build_script, "--input", input_path, "--output", output_path, "--min-samples-per-group", "3"),
    env = c("OMP_NUM_THREADS=1", "R_MAX_VSIZE=100Gb"),
    stdout = "",
    stderr = ""
  )

  if (!identical(status, 0L)) {
    stop("Diffuse-vs-Intestinal DEG cache build failed for: ", input_path)
  }
}

message("")
message("Finished checking/building Diffuse-vs-Intestinal DEG caches.")
