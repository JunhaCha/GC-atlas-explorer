#!/usr/bin/env Rscript

atlas_paths <- c(
  "/Users/junhacha/Documents/Playground/seurat_merged_TME_malignant_final_umap_app_slim.rds",
  "/Users/junhacha/Documents/Playground/seurat_epithelial_normal_final_final_app_slim.rds",
  "/Users/junhacha/Documents/Playground/seurat_cancercells_final_app_slim.rds",
  "/Users/junhacha/Documents/Playground/seurat_Stromal_final_app_slim.rds",
  "/Users/junhacha/Documents/Playground/seurat_CD8T_final2_app_slim.rds",
  "/Users/junhacha/Documents/Playground/seurat_CD4T_final2_app_slim.rds",
  "/Users/junhacha/Documents/Playground/seurat_B_final2_app_slim.rds",
  "/Users/junhacha/Documents/Playground/seurat_Mye_final2_app_slim.rds"
)

build_script <- "/Users/junhacha/Documents/Playground/seurat-shiny-explorer/build_average_cache_for_app.R"
force_rebuild <- "--force" %in% commandArgs(trailingOnly = TRUE)

for (input_path in atlas_paths) {
  if (!file.exists(input_path)) {
    message("Skipping missing slim object: ", input_path)
    next
  }

  output_path <- sub("\\.rds$", "_avg_cache.rds", input_path, ignore.case = TRUE)
  if (file.exists(output_path) && !force_rebuild) {
    message("Average cache already exists, skipping: ", output_path)
    next
  }

  message("")
  message("==== Building Average Cache ====")
  message("Input:  ", input_path)
  message("Output: ", output_path)

  status <- system2(
    command = Sys.which("Rscript"),
    args = c(build_script, "--input", input_path, "--output", output_path),
    env = c("OMP_NUM_THREADS=1", "R_MAX_VSIZE=100Gb"),
    stdout = "",
    stderr = ""
  )

  if (!identical(status, 0L)) {
    stop("Average-cache build failed for: ", input_path)
  }
}

message("")
message("Finished checking/building average-expression caches.")
