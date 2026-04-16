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
  "seurat_merged_TME_malignant_final_umap.rds",
  "seurat_epithelial_normal_final_final.rds",
  "seurat_cancercells_final.rds",
  "seurat_Stromal_final.rds",
  "seurat_CD8T_final2.rds",
  "seurat_CD4T_final2.rds",
  "seurat_B_final2.rds",
  "seurat_Mye_final2.rds"
)
atlas_paths <- file.path(app_data_dir(), atlas_filenames)

slim_script <- file.path(repo_dir, "slim_seurat_for_app.R")

for (input_path in atlas_paths) {
  if (!file.exists(input_path)) {
    message("Skipping missing file: ", input_path)
    next
  }

  output_path <- sub("\\.rds$", "_app_slim.rds", input_path, ignore.case = TRUE)
  if (file.exists(output_path)) {
    message("Slim file already exists, skipping: ", output_path)
    next
  }

  cmd <- c(
    slim_script,
    "--input", input_path,
    "--output", output_path
  )

  message("")
  message("==== Slimming ====")
  message("Input:  ", input_path)
  message("Output: ", output_path)

  status <- system2(
    command = Sys.which("Rscript"),
    args = cmd,
    env = c("OMP_NUM_THREADS=1", "R_MAX_VSIZE=100Gb"),
    stdout = "",
    stderr = ""
  )

  if (!identical(status, 0L)) {
    stop("Slimming failed for: ", input_path)
  }
}

message("")
message("Finished checking/building slim Seurat files.")
