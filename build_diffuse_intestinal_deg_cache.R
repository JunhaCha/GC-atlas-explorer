#!/usr/bin/env Rscript

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  parsed <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected argument: %s", key))
    }
    key <- sub("^--", "", key)
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      parsed[[key]] <- TRUE
      i <- i + 1
    } else {
      parsed[[key]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  parsed
}

get_arg <- function(args, name, default = NULL, required = FALSE) {
  value <- args[[name]]
  if (is.null(value)) {
    if (required) {
      stop(sprintf("Missing required argument --%s", name))
    }
    return(default)
  }
  value
}

args <- parse_args()
input_path <- get_arg(args, "input", required = TRUE)
output_path <- get_arg(args, "output", required = TRUE)
group_col <- get_arg(args, "group-col", default = "final_group")
cluster_col <- get_arg(args, "cluster-col", default = "final_celltype")
min_samples_per_group <- as.integer(get_arg(args, "min-samples-per-group", default = 3))

script_path <- normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]),
  winslash = "/",
  mustWork = FALSE
)
repo_dir <- dirname(script_path)
source(file.path(repo_dir, "R", "helpers.R"))

message("Reading: ", input_path)
avg_cache_path <- sub("\\.rds$", "_avg_cache.rds", input_path, ignore.case = TRUE)
if (!file.exists(avg_cache_path)) {
  stop("Average-expression cache not found: ", avg_cache_path)
}
avg_cache <- readRDS(avg_cache_path)

cache <- build_diffuse_intestinal_deg_cache(
  avg_cache = avg_cache,
  group_col = group_col,
  cluster_col = cluster_col,
  min_samples_per_group = min_samples_per_group
)

message("Writing: ", output_path)
saveRDS(cache, output_path, compress = TRUE)
message(sprintf("Cache rows: %d", nrow(cache$deg_table)))
