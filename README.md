# Gastric Cancer single-cell Explorer

Shiny app for exploring gastric cancer single-cell atlas objects plus lightweight Xenium spatial exports, precomputed marker tables, average-expression caches, and malignant-vs-normal quadrant DEG plots.

## Runtime layout

The app is now runtime-configurable and does not require Mac-specific paths.

- Code directory: this repository
- Data directory: controlled by `GC_APP_DATA_DIR`
- Default local behavior: if `GC_APP_DATA_DIR` is unset, the app looks for data in `./data` and then in the parent directory of the app

## Required runtime files

Place the files from data directory mounted into the container (GCshinyapp_data folder):

- `seurat_merged_TME_malignant_final_umap_app_slim.rds`
- `seurat_epithelial_normal_final_final_app_slim.rds`
- `seurat_cancercells_final_app_slim.rds`
- `seurat_Stromal_final_app_slim.rds`
- `seurat_CD8T_final2_app_slim.rds`
- `seurat_CD4T_final2_app_slim.rds`
- `seurat_B_final2_app_slim.rds`
- `seurat_Mye_final2_app_slim.rds`
- `*_avg_cache.rds` files for the objects you want to support
- `*_quadrant_cache.rds` for quadrant DEG plots
- `dataset_group_colors.rds`
- `path_group_colors.rds`
- `subcluster_ct_colors_combined.rds`
- `markers_RNA_merged_AllCelltypes_final.xlsx`
- marker/DE `.xlsx` files for the other atlas objects you want exposed
- `atlas_app_inputs/export_manifest.tsv`
- `atlas_app_inputs/<sample_id>/cells.rds`
- `atlas_app_inputs/<sample_id>/image_lowres.png`
- `atlas_app_inputs/<sample_id>/image_meta.rds`
- `atlas_app_inputs/<sample_id>/neighborhood_summary.rds`
- `atlas_app_inputs/<sample_id>/neighborhood_composition.rds`
- `atlas_app_inputs/<sample_id>/<sample_id>_counts_sce.rds` if you want to rebuild Xenium marker assets
- `atlas_app_inputs/<sample_id>/marker_expr/` if you want the Xenium Marker Overlay tab available immediately

Example shared host path used for this project:

- `/mnt/ix1/Projects_lite/junhacha/GCshinyapp_data`

## Local run

```bash
cd /path/to/GC-atlas-explorer
GC_APP_DATA_DIR=/path/to/data R -e "shiny::runApp('.', launch.browser = TRUE)"
```

For the shared project storage, that becomes:

```bash
GC_APP_DATA_DIR=/mnt/ix1/Projects_lite/junhacha/GCshinyapp_data R -e "shiny::runApp('.', launch.browser = TRUE)"
```

## Docker run

Build:

```bash
docker build -t gc-singlecell-explorer .
```

Run:

```bash
docker run --rm -p 3838:3838 \
  -e GC_APP_DATA_DIR=/data \
  -v /srv/gc-explorer-data:/data:ro \
  gc-singlecell-explorer
```

Then open:

- `http://<server-host>:3838`

## Docker Compose

Example file is included at [docker-compose.yml](/Users/junhacha/Documents/Playground/GC-atlas-explorer/docker-compose.yml).

Set the host data path in `.env`:

```bash
cp .env.example .env
```

Then edit:

```bash
GC_APP_HOST_DATA_DIR=/srv/gc-explorer-data
```

Start:

```bash
docker compose up -d --build
```

## Admin handoff
1. This repository
2. A host data directory containing the required `.rds`, cache, and `.xlsx` files
3. One of:
   - `docker build` + `docker run`
   - `docker compose up -d --build`

## Updating the app later

You can keep extending the app after deployment. Update the code in this repo, rebuild the container, and restart it. The shared data directory can stay mounted in the same place.
