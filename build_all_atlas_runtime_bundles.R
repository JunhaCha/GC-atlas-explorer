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
  "seurat_merged_TME_malignant_final_umap_app_slim.rds",
  "seurat_epithelial_normal_final_final_app_slim.rds",
  "seurat_cancercells_final_app_slim.rds",
  "seurat_Stromal_final_app_slim.rds",
  "seurat_CD8T_final2_app_slim.rds",
  "seurat_CD4T_final2_app_slim.rds",
  "seurat_B_final2_app_slim.rds",
  "seurat_Mye_final2_app_slim.rds"
)

build_script <- file.path(repo_dir, "build_atlas_runtime_bundle.R")

for (input_path in file.path(app_data_dir(), atlas_filenames)) {
  if (!file.exists(input_path)) {
    message("Skipping missing slim object: ", input_path)
    next
  }
  message("")
  message("==== Building Atlas Runtime Bundle ====")
  message("Input: ", input_path)
    status <- system2(
      command = "/Library/Frameworks/R.framework/Resources/bin/Rscript",
      args = c(build_script, "--input", input_path),
      stdout = "",
      stderr = ""
    )
    if (!identical(status, 0L)) {
        message("Skipping bundle build after failure: ", input_path)
        next
    }
}

message("")
message("Finished building atlas runtime bundles.")
