#!/usr/bin/env Rscript

atlas_paths <- c(
  "/Users/junhacha/Documents/Playground/seurat_merged_TME_malignant_final_umap.rds",
  "/Users/junhacha/Documents/Playground/seurat_epithelial_normal_final_final.rds",
  "/Users/junhacha/Documents/Playground/seurat_cancercells_final.rds",
  "/Users/junhacha/Documents/Playground/seurat_Stromal_final.rds",
  "/Users/junhacha/Documents/Playground/seurat_CD8T_final2.rds",
  "/Users/junhacha/Documents/Playground/seurat_CD4T_final2.rds",
  "/Users/junhacha/Documents/Playground/seurat_B_final2.rds",
  "/Users/junhacha/Documents/Playground/seurat_Mye_final2.rds"
)

slim_script <- "/Users/junhacha/Documents/Playground/seurat-shiny-explorer/slim_seurat_for_app.R"

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
