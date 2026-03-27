FROM rocker/shiny:4.4.1

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libglpk-dev \
    libhdf5-dev \
    make \
    g++ \
    git \
    curl \
  && rm -rf /var/lib/apt/lists/*

RUN install2.r --error --skipinstalled \
    shiny \
    bslib \
    dplyr \
    ggplot2 \
    viridisLite \
    pheatmap \
    readxl \
    ggrepel \
    tibble \
    Seurat

WORKDIR /opt/gc-singlecell-explorer
COPY . /opt/gc-singlecell-explorer

ENV GC_APP_DATA_DIR=/data
EXPOSE 3838

CMD ["R", "-e", "options(shiny.host='0.0.0.0', shiny.port=3838); shiny::runApp('/opt/gc-singlecell-explorer', launch.browser = FALSE)"]
