# Gastric Cancer single-cell Explorer

Shiny app for exploring gastric cancer single-cell Seurat objects, precomputed marker tables, average-expression caches, and malignant-vs-normal quadrant DEG plots.

## Runtime layout

The app is now runtime-configurable and does not require Mac-specific paths.

- Code directory: this repository
- Data directory: controlled by `GC_APP_DATA_DIR`
- Default local behavior: if `GC_APP_DATA_DIR` is unset, the app looks for data in `./data` and then in the parent directory of the app

## Required runtime files

Place the files from data directory mounted into the container (GCshinyapp_data folder):



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

## to admin (Lucas)
1. This repository
2. A host data directory containing the required `.rds`, cache, and `.xlsx` files
3. One of:
   - `docker build` + `docker run`
   - `docker compose up -d --build`

## Updating the app later

will update later to include Xenium.

