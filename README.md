# Gastric Cancer single-cell Explorer

Shiny app for exploring gastric cancer single-cell Seurat objects, precomputed marker tables, average-expression caches, and malignant-vs-normal quadrant DEG plots.

## Runtime layout

The app is now runtime-configurable and does not require Mac-specific paths.

- Code directory: this repository
- Data directory: controlled by `GC_APP_DATA_DIR`
- Default local behavior: if `GC_APP_DATA_DIR` is unset, the app looks for data in `./data` and then in the parent directory of the app

## Required runtime files

Place these in the data directory mounted into the container or available on the host:

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

## Local run

```bash
cd /path/to/seurat-shiny-explorer
GC_APP_DATA_DIR=/path/to/data R -e "shiny::runApp('.', launch.browser = TRUE)"
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

Example file is included at [docker-compose.yml](/Users/junhacha/Documents/Playground/seurat-shiny-explorer/docker-compose.yml).

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

The administrator only needs:

1. This repository
2. A host data directory containing the required `.rds`, cache, and `.xlsx` files
3. One of:
   - `docker build` + `docker run`
   - `docker compose up -d --build`

## Updating the app later

Yes, you can absolutely keep modifying the app after this.

- Change code in this repository
- Rebuild the container image
- Restart the container

The data can stay mounted in the same host directory, so code updates do not require rebundling the large Seurat objects every time.
