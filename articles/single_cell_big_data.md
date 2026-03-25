# Scaling to millions of cells with bixverse

## Intro

This vignette demonstrates how to run a single cell analysis on a data
set containing nearly 3 million cells on a local machine - in this case
a MacBook Pro M1 Max - without exhausting memory. If you have not read
the [design
choices](https://gregorlueg.github.io/bixverse/docs/articles/design_single_cell.html)
and the [introductory
vignette](https://gregorlueg.github.io/bixverse/docs/articles/single_cell_intro.html),
please do so first. The [PBMC3k
walkthrough](https://gregorlueg.github.io/bixverse/docs/articles/pbmc3k.html)
covers every step in detail on a small data set; here we focus on what
changes when the cell count moves from thousands to millions.

The data set used is one of the h5ad of the [100 million cell
Tahoe](https://www.biorxiv.org/content/10.1101/2025.02.20.639398v3) data
set from Parse Biosciences. The original file has 5.5m cells; after
quality control roughly 2.9 million cells are retained and read in.

> **Note**
>
> None of the code chunks below are evaluated. This vignette is a
> reference for how one *would* run such an analysis, not a live
> rendering - the data set is far too large for CI/CD runners and
> Microsoft would not be happy.

If you want to reproduce this, download one of the h5ad files or data
you always wanted to analyse.

``` r
library(bixverse)
library(ggplot2)
library(data.table)
```

## Loading the data

For large `.h5ad` files the key entry point is `stream_h5ad`, which
streams the sparse matrix into the Rust binary files and the internal
DuckDB. Once the data has been streamed once, subsequent sessions can
reconnect via `load_existing`.

``` r
h5_path <- path.expand(
  "PATH TO H5AD/plate1_filt_Vevo_Tahoe100M_WServicesFrom_ParseGigalab.h5ad"
)

dir_files <- "PATH TO WHERE TO SAVE THE FILES"

# Quick sanity check on dimensions before committing to the stream
meta <- get_h5ad_dimensions(h5_path)
meta

sc_object <- SingleCells(
  dir_data = dir_files
)

sc_object <- stream_h5ad(object = sc_object, h5_path = h5_path)
```

This step is slow. Ca. 5 to 10 minutes. If you have already loaded in
the data and exist on-disk from a previous session, skip the streaming
step entirely and reconnect:

``` r
sc_object <- SingleCells(
  dir_data = dir_files
)

sc_object <- load_existing(object = sc_object)
```

## Quality control

### Gene set proportions

The same gene set proportion logic from the PBMC3k vignette applies, but
the `streaming` argument **must** be set to `TRUE` so that proportions
are computed in batches rather than pulling the entire count matrix into
Rust… The memory pressure you should observe should remain minimal.

``` r
var <- get_sc_var(object = sc_object)

gs_of_interest <- list(
  MT = var[grepl("^MT-", gene_id), gene_id],
  Ribo = var[grepl("^RPS|^RPL", gene_id), gene_id]
)

sc_object <- gene_set_proportions_sc(
  object = sc_object,
  gene_set_list = gs_of_interest,
  streaming = TRUE,
  .verbose = TRUE
)

sc_object[[1:5L]]
```

### MAD outlier detection

Outlier detection works identically to the small-data case. The QC
metrics table is modest in size (one row per cell, a handful of columns)
and fits comfortably in memory even at millions of cells.

``` r
qc_df <- sc_object[[c("cell_id", "lib_size", "nnz", "MT")]]

metrics <- list(
  log10_lib_size = log10(qc_df$lib_size),
  log10_nnz = log10(qc_df$nnz),
  MT = qc_df$MT
)

directions <- c(
  log10_lib_size = "twosided",
  log10_nnz = "twosided",
  MT = "above"
)

qc <- run_cell_qc(metrics, directions, threshold = 3)

# not needed anymore
rm(qc_df)

qc
```

### Filtering cells

``` r
sc_object[["outlier"]] <- qc$combined

cells_to_keep <- qc_df[!qc$combined, cell_id]

sc_object <- set_cells_to_keep(sc_object, cells_to_keep)

sc_object
```

## Feature selection and PCA

HVG selection at this scale again requires `streaming = TRUE`. PCA uses
a sparse SVD, which is automatically selected when cell counts exceed
500k but can be set explicitly.

``` r
sc_object <- find_hvg_sc(
  object = sc_object,
  streaming = TRUE,
  .verbose = TRUE
)

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = 30L,
  sparse_svd = TRUE
)
```

## Nearest neighbours (GPU-accelerated)

At nearly 3 million cells, nearest neighbour search is the main
computational bottleneck. The companion package
[bixverse.gpu](https://github.com/GregorLueg/bixverse.gpu) exposes a
CAGRA-style GPU-accelerated approximate nearest neighbour search that
makes this quite fast on a laptop with a discrete GPU - thanks to CubeCL
and WGPU as long as you can fit the data into VRAM (should be feasible
with 8 GBs - I think?) this should run quite fast. CPU-side alternatives
such as `nndescent` or `HNSW` also work but will be considerably slower
at this scale. For reference:

- GPU version on M1 Max with 64 GB takes ~1 minute.
- HNSW somewhere between 3 to 4 minutes.

``` r
library(bixverse.gpu)

sc_object <- find_neighbours_cagra_sc(
  object = sc_object,
  cagra_params = params_sc_cagra(k_query = 15L, ann_dist = "cosine"),
  extract_knn = FALSE,
  .verbose = TRUE
)
```

## Visualisation

Running UMAP on the shared nearest neighbour graph built from ~2.9
million cells and 65 million+ edges takes time but does complete on a
single machine (optimisation with parallel Adam optimisation done in \<3
minutes for me).

``` r
sc_object <- umap_sc(object = sc_object)
```

## Differential gene expression

With metadata stored in the obs table, subsetting cells for a pairwise
comparison is straightforward. Here we compare two drug treatment
conditions. Again, we leverage the streaming here… Instead of holding
the full data set in memory, we just load in the two groups (each
roughly 30k cells), do a fast ranking in Rust and do the DGE. Easy
peezy, done in ca. 3 seconds.

``` r
meta_data <- sc_object[[c("cell_id", "cell_name", "drug")]]

cells_grp_a <- meta_data[drug == "AT7519", cell_id]
cells_grp_b <- meta_data[drug == "palbociclib", cell_id]

dge_results <- find_markers_sc(
  object = sc_object,
  cells_1 = cells_grp_a,
  cells_2 = cells_grp_b
)

setorder(dge_results, fdr)

head(dge_results, 25L)
```

## Key differences from the small-data workflow

The API surface is intentionally identical to the PBMC3k walkthrough.
The meaningful differences when working at the million-cell scale are:

1.  **Streaming I/O.** Use `stream_h5ad` instead of `load_h5ad`, and set
    `streaming = TRUE` on `gene_set_proportions_sc` and `find_hvg_sc`.
    This is VERY important.
2.  **On-disk persistence.** `load_existing` lets you reconnect to a
    previously streamed data set without re-reading the `.h5ad`. (Btw,
    works for any data.)
3.  **Sparse SVD.** Automatically enabled above 500k cells; can be
    forced explicitly via `sparse_svd = TRUE`.
4.  **GPU-accelerated kNN.** `bixverse.gpu` provides CAGRA-based nearest
    neighbour search that scales comfortably to millions of cells and
    makes it run in a minute or so. The CPU versions are still highly
    optimised with SIMD and specifically designed memory layouts to
    avoid making your memory cry.
