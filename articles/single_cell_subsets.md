# Sub-clustering and per-group pipelines

## Intro

Your first pass of clustering gives you lineages: T cells, B cells,
monocytes. What it does not give you is the structure inside those
lineages. The highly variable genes and the principal components were
chosen to separate a monocyte from a B cell, and on that axis every T
cell looks the same. If you crank up the Leiden resolution on the full
object you mostly get the lineages chopped into arbitrary pieces, not
naive versus memory CD4.

The fix is to throw the global feature selection away and redo it on the
cells you actually care about. That is what `SingleCellsSubset` is for.
It gives you a view onto a subset of a
[`SingleCells`](https://gregorlueg.github.io/bixverse/reference/SingleCells.html)
object that shares the on-disk counts with the parent, but carries its
own empty cache, so HVG, PCA, neighbours and clusters all get recomputed
from scratch on the subset alone.

This vignette assumes you have read the [design
choices](https://gregorlueg.github.io/bixverse/articles/design_single_cell.html)
and the [PBMC
walkthrough](https://gregorlueg.github.io/bixverse/articles/pbmc_single_cell.html).
The standard chain is not re-explained here.

``` r

library(bixverse)
library(ggplot2)
library(data.table)
#> 
#> Attaching package: 'data.table'
#> The following object is masked from 'package:base':
#> 
#>     %notin%
library(bixverse.plots)
library(magrittr)
```

## A coarse first pass

Same PBMC3k data set as the main walkthrough, run through QC, HVG, PCA,
neighbours and Leiden without commentary.

``` r

pbmc3k_path <- download_pbmc3k()

tempdir_pbmc <- tempdir()

sc_object <- SingleCells(dir_data = tempdir_pbmc)

sc_object <- load_mtx(
  object = sc_object,
  sc_mtx_io_param = get_cell_ranger_params(pbmc3k_path),
  mtx_streaming = FALSE,
  .verbose = FALSE
)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

setnames_sc(
  object = sc_object,
  table = "var",
  old = "column1",
  new = "gene_symbol"
)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

var <- get_sc_var(sc_object)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

ensembl_to_symbol <- setNames(var$gene_symbol, var$gene_id)
symbol_to_ensembl <- setNames(var$gene_id, var$gene_symbol)
```

Mitochondrial proportions, MAD outliers, cells to keep.

``` r

gs_of_interest <- list(
  MT = var[grepl("^MT-", gene_symbol), gene_id]
)

sc_object <- gene_set_proportions_sc(
  sc_object,
  gs_of_interest,
  streaming = FALSE,
  .verbose = FALSE
)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

qc_df <- sc_object[[c("cell_id", "lib_size", "nnz", "MT")]]
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

qc <- run_cell_qc(
  metrics = list(
    log10_lib_size = log10(qc_df$lib_size),
    log10_nnz = log10(qc_df$nnz),
    MT = qc_df$MT
  ),
  cells_to_keep = get_cells_to_keep(sc_object),
  directions = c(
    log10_lib_size = "twosided",
    log10_nnz = "twosided",
    MT = "above"
  ),
  threshold = 3
)

sc_object <- set_cells_to_keep(sc_object, qc_df[!qc$combined, cell_id])
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
```

Now the standard chain. The data is tiny, so exhaustive kNN beats
building an index.

``` r

sc_object <- find_hvg_sc(sc_object, hvg_no = 2000L, .verbose = FALSE)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = 30L,
  sparse_svd = TRUE,
  .verbose = FALSE
)

sc_object <- find_neighbours_sc(
  object = sc_object,
  neighbours_params = params_sc_neighbours(
    knn = list(knn_method = "exhaustive")
  ),
  .verbose = FALSE
)

sc_object <- find_clusters_sc(sc_object, res = 1, name = "leiden_global")
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

sc_object <- umap_sc(sc_object, .verbose = FALSE)
```

For the annotation we deliberately stay coarse. No CD4 versus CD8 split,
no naive versus memory. Lineage only, because the fine structure is
exactly what we want the subset to recover.

``` r

coarse_markers <- c(
  CD3D = "T cells",
  CD3E = "T cells",
  CD3G = "T cells",
  TRAC = "T cells",
  MS4A1 = "B cells",
  CD79A = "B cells",
  CD19 = "B cells",
  CD14 = "Monocytes",
  LYZ = "Monocytes",
  S100A8 = "Monocytes",
  FCGR3A = "Monocytes",
  GNLY = "NK",
  NKG7 = "NK",
  KLRD1 = "NK",
  FCER1A = "DC",
  CD1C = "DC",
  LILRA4 = "DC",
  PPBP = "Platelet",
  PF4 = "Platelet"
)

coarse_dt <- stack(coarse_markers) %>%
  as.data.table() %>%
  setnames(., c("values", "ind"), c("cell_type", "gene_symbol")) %>%
  .[, gene_symbol := as.character(gene_symbol)] %>%
  .[, gene_id := var$gene_id[match(gene_symbol, var$gene_symbol)]] %>%
  .[!is.na(gene_id), ]

sctype_scores <- calc_sc_type_scores(
  object = sc_object,
  cell_marker_list = prepare_cell_markers(sc_object, coarse_dt)
)

cell_type_anno <- score_clusters(
  sctype_scores,
  sc_object[[]][["leiden_global"]]
)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

obs <- get_sc_obs(sc_object, filtered = TRUE)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

sc_object[["cell_type"]] <- cell_type_anno$cell_type[match(
  obs$leiden_global,
  cell_type_anno$cluster_id
)]
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

get_sc_obs(sc_object, filtered = TRUE)[, .N, by = cell_type][order(-N)]
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#>    cell_type     N
#>       <char> <int>
#> 1:   T cells  1043
#> 2: Monocytes   459
#> 3:        NK   382
#> 4:   B cells   279
```

``` r

embedding_plot_sc(
  sc_object,
  embedding = "umap",
  colour_by = "cell_type",
  label_by = "cell_type",
  discrete = TRUE
) +
  labs(colour = "Cell type")
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
```

![](single_cell_subsets_files/figure-html/global%20umap-1.png)

## Taking a subset

The constructor takes the parent, a column in obs, and the level you
want.

``` r

t_cells <- SingleCellsSubset(
  sc_object = sc_object,
  grouping_column = "cell_type",
  group = "T cells"
)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

t_cells
#> Single cell experiment (subset).
#>   No cells: 1043
#>   No genes: 11139
#>   Group: cell_type = T cells
#>   HVG calculated: FALSE
#>   PCA calculated: FALSE
#>   Other embeddings: none
#>   KNN generated: FALSE
#>   SNN generated: FALSE
```

The print method tells you the important bit: HVG, PCA, kNN and sNN are
all `FALSE`. Nothing was inherited. The subset shares the Rust count
connection with the parent, so no counts were copied and no files were
written, but the `ScCache` is brand new and the parent’s HVG selection
was dropped on purpose.

Everything else behaves like the parent object. The obs table is
filtered to the group, the var table is untouched, and the usual getters
work.

``` r

dim(t_cells)
#> [1]  1043 11139

head(t_cells)[, .(cell_id, cell_idx, lib_size, cell_type)]
#>             cell_id cell_idx lib_size cell_type
#>              <char>    <num>    <num>    <char>
#> 1: AAACATACAACCAC-1        1     2418   T cells
#> 2: AAACATTGATCAGC-1        3     3144   T cells
#> 3: AAACGCACTGGTAC-1        6     2154   T cells
#> 4: AAACTTGATCCAGA-1       12     2382   T cells
#> 5: AAAGAGACGAGATA-1       13     2403   T cells
#> 6: AAAGCCTGTATGCG-1       20     2922   T cells
```

### A note on indices

This is the one place where you can shoot yourself in the foot, so it is
worth being explicit. The subset keeps two index spaces around:

- **Subset positions** are 1-indexed rows of the subset’s own obs table.
  This is what `[` and `[[` take.
- **Original positions** are the parent’s index space.
  `subset_to_original` and the `cell_idx` obs column both hold these,
  1-indexed. Everything that reaches Rust gets translated into this
  space, 0-indexed, before the call.

[`get_cell_indices()`](https://gregorlueg.github.io/bixverse/reference/get_cell_indices.md)
lets you ask for either.

``` r

first_three <- get_cell_names(t_cells)[1:3]

# 1-indexed position within the subset
get_cell_indices(t_cells, cell_ids = first_three, rust_index = FALSE)
#> [1] 1 2 3

# 0-indexed position in the parent, which is what Rust sees
get_cell_indices(t_cells, cell_ids = first_three, rust_index = TRUE)
#> [1] 0 2 5
```

The practical upshot: you never translate manually. Hand cell IDs or
subset positions to the getters and let the class deal with it. The
`cell_idx` column is there so you can join results back onto the parent
later, which is exactly what we do further down.

## Re-running the chain

Same functions as on the parent. They dispatch on the subset.

``` r

t_cells <- find_hvg_sc(t_cells, hvg_no = 2000L, .verbose = FALSE)

t_cells <- calculate_pca_sc(
  object = t_cells,
  no_pcs = 20L,
  sparse_svd = TRUE,
  .verbose = FALSE
)

t_cells <- find_neighbours_sc(
  object = t_cells,
  neighbours_params = params_sc_neighbours(
    knn = list(knn_method = "exhaustive")
  ),
  .verbose = FALSE
)

t_cells <- find_clusters_sc(t_cells, res = 0.5, name = "t_subcluster")

t_cells <- umap_sc(t_cells, .verbose = FALSE)

t_cells
#> Single cell experiment (subset).
#>   No cells: 1043
#>   No genes: 11139
#>   Group: cell_type = T cells
#>   HVG calculated: TRUE
#>   PCA calculated: TRUE
#>   Other embeddings: umap
#>   KNN generated: TRUE
#>   SNN generated: TRUE
```

Was any of that worth it? Check how much the feature selection actually
moved.

``` r

global_hvg <- get_hvg(sc_object)
subset_hvg <- get_hvg(t_cells)

length(intersect(global_hvg, subset_hvg))
#> [1] 1021
```

Around a thousand of the 2000 genes are different, and the ones that
dropped out are the monocyte and B cell markers that were carrying the
global PCA. The subset PCA is built on genes that vary *within* T cells,
which is the whole point.

``` r

embedding_plot_sc(
  t_cells,
  embedding = "umap",
  colour_by = "t_subcluster",
  label_by = "t_subcluster",
  discrete = TRUE
) +
  labs(colour = "Sub-cluster")
```

![](single_cell_subsets_files/figure-html/subset%20umap-1.png)

## Markers for the sub-clusters

Differential expression dispatches on the subset, so the comparison is
automatically confined to it. Each sub-cluster is tested against the
other T sub-clusters, not against every monocyte in the data set.

``` r

sub_markers <- find_all_markers_sc(
  object = t_cells,
  column_of_interest = "t_subcluster",
  .verbose = FALSE
)

sub_markers[, gene_symbol := ensembl_to_symbol[gene_id]]

sub_markers[fdr <= 0.05][order(grp, -lfc)][,
  head(.SD, 5),
  by = grp,
  .SDcols = c("gene_symbol", "lfc", "fdr")
]
#>       grp gene_symbol       lfc          fdr
#>     <int>      <char>     <num>        <num>
#>  1:     0      S100A4 1.4256958 7.233787e-73
#>  2:     0        KLF6 0.8074021 7.268987e-25
#>  3:     0       ANXA1 0.7887582 2.502905e-24
#>  4:     0        IL32 0.7337589 7.998436e-26
#>  5:     0       CLIC1 0.7258839 1.613255e-23
#>  6:     1        CCR7 0.5556238 3.884105e-15
#>  7:     1     C6orf48 0.4032111 8.996014e-10
#>  8:     1        CD8B 0.3322895 2.215998e-03
#>  9:     1        SELL 0.3156452 1.354104e-04
#> 10:     1        BTG1 0.3132365 3.699444e-10
```

S100A4, ANXA1 and IL32 on one side, CCR7, SELL and CD8B on the other.
That is the memory versus naive split, and the global clustering had no
way of seeing it. Cluster 2 is two cells, which is Leiden being Leiden
at this resolution.

The expression plots work on the subset too, so you can go straight to a
dot plot without touching the parent.

``` r

t_features <- setdiff(
  symbol_to_ensembl[c(
    "CCR7",
    "LEF1",
    "SELL",
    "IL7R",
    "S100A4",
    "IL32",
    "CD8A",
    "CD8B",
    "GZMK",
    "NKG7"
  )],
  NA
)

dot_plot_sc(
  object = t_cells,
  features = t_features,
  feature_labels = ensembl_to_symbol[t_features],
  grouping_variable = "t_subcluster",
  scale_exp = TRUE,
  cluster_groups = TRUE
)
```

![](single_cell_subsets_files/figure-html/subcluster%20dotplot-1.png)

## Getting the labels back to the parent

So far everything lives in the subset, in memory. To make the labels
permanent, or to plot them on the parent’s global UMAP, push them into
the parent’s DuckDB.
[`merge_subset_obs()`](https://gregorlueg.github.io/bixverse/reference/merge_subset_obs.md)
does the join on `cell_idx` for you.

``` r

sc_object <- merge_subset_obs(
  object = sc_object,
  subsets = t_cells,
  cols = "t_subcluster",
  prefix_values = TRUE
)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> Merged 1 column(s) for 1043 cells into obs (1657 cells left as NA).

get_sc_obs(sc_object, filtered = TRUE)[, .N, by = t_subcluster][order(
  t_subcluster
)]
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#>    t_subcluster     N
#>          <char> <int>
#> 1:    T cells_0   562
#> 2:    T cells_1   481
#> 3:         <NA>  1120
```

`prefix_values = TRUE` stamps the subset’s group onto every value, so
you get `T cells_0` rather than a bare `0`. Drop it if you want the raw
labels. Leave `cols` out and everything the subset gained that the
parent does not have yet gets merged, which after a pipeline run is
exactly the new columns.

The join is a left join, so every cell that was not a T cell gets `NA`.
That is the behaviour you want, and the downstream methods respect it.

Results that already carry their own `cell_idx` skip the helper
entirely.
[`fast_cluster_sc()`](https://gregorlueg.github.io/bixverse/reference/fast_cluster_sc.md)
returns memberships keyed in parent space, so
[`add_sc_new_obs()`](https://gregorlueg.github.io/bixverse/reference/add_sc_new_obs.md)
takes them directly.

``` r

fast_sub <- fast_cluster_sc(
  object = t_cells,
  resolutions = c(1, 0.5),
  .verbose = FALSE
)

sc_object <- add_sc_new_obs(sc_object, get_data(fast_sub))
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

head(sc_object)[, .(cell_id, cell_type, t_subcluster, res_1, res_0.5)]
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#>             cell_id cell_type t_subcluster res_1 res_0.5
#>              <char>    <char>       <char> <int>   <int>
#> 1: AAACATACAACCAC-1   T cells    T cells_0     0       0
#> 2: AAACATTGATCAGC-1   T cells    T cells_0     0       0
#> 3: AAACCGTGCTTCCG-1 Monocytes         <NA>    NA      NA
#> 4: AAACGCACTGGTAC-1   T cells    T cells_0     0       0
#> 5: AAACGCTGGTTCTT-1        NK         <NA>    NA      NA
#> 6: AAACGCTGTAGCCA-1        NK         <NA>    NA      NA
```

## Pipelines

Once you are doing this for more than one cell type, writing out HVG,
PCA, neighbours and clusters by hand five times gets old.
[`sc_pipeline()`](https://gregorlueg.github.io/bixverse/reference/sc_pipeline.md)
lets you declare the chain once.

``` r

subcluster_pipeline <- sc_pipeline() %>>%
  step_hvg_sc(hvg_no = 2000L, .verbose = FALSE) %>>%
  step_pca_sc(no_pcs = 20L, sparse_svd = TRUE, .verbose = FALSE) %>>%
  step_neighbours_sc(
    neighbours_params = params_sc_neighbours(
      knn = list(knn_method = "exhaustive")
    ),
    .verbose = FALSE
  ) %>>%
  step_clusters_sc(res = 0.5, name = "subcluster")

subcluster_pipeline
#> <ScPipeline> 4 steps
#>   1. hvg         hvg_no = 2000L, hvg_params = <list>, streaming = NULL, .verbose = FALSE
#>   2. pca         no_pcs = 20L, pca_params = <list>, sparse_svd = TRUE, hvg = NULL, seed = 42L, .verbose = FALSE
#>   3. neighbours  embd_to_use = "pca", no_embd_to_use = NULL, modality = c("rna", "adt"), neighbours_params = <list>, seed = 42L, .verbose = FALSE
#>   4. clusters    cluster_algorithm = c("leiden", "louvain"), res = 0.5, name = "subcluster", modality = c("rna", "adt", "wnn"), seed = 42L
```

Pipelines are inert. Nothing runs until you apply one, and the same
pipeline works on a `SingleCells` or a `SingleCellsSubset` because
dispatch happens inside each step.

``` r

b_cells <- SingleCellsSubset(
  sc_object = sc_object,
  grouping_column = "cell_type",
  group = "B cells"
)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpSlHbTJ/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> ℹ This message has been shown 60 times and will not be shown again this session.

b_cells <- apply_pipeline(subcluster_pipeline, b_cells)

b_cells
#> Single cell experiment (subset).
#>   No cells: 279
#>   No genes: 11139
#>   Group: cell_type = B cells
#>   HVG calculated: TRUE
#>   PCA calculated: TRUE
#>   Other embeddings: none
#>   KNN generated: TRUE
#>   SNN generated: TRUE
```

### One pipeline, every group

[`apply_pipeline_per_group()`](https://gregorlueg.github.io/bixverse/reference/apply_pipeline_per_group.md)
builds the subset for you and returns a named list.

Sizing matters here. PCA with 20 components needs more than 20 cells, so
a rare population that the annotation happened to call will fall over
and take the whole run with it. Restrict the groups explicitly rather
than hoping for the best. On this data set the coarse annotation only
returned four lineages and all of them clear the bar, so nothing gets
dropped, but the guard costs nothing.

``` r

cell_counts <- get_sc_obs(sc_object, filtered = TRUE)[, .N, by = cell_type][
  order(-N)
]

groups_to_run <- cell_counts[N >= 100, cell_type]

groups_to_run
#> [1] "T cells"   "Monocytes" "NK"        "B cells"
```

``` r

per_group <- apply_pipeline_per_group(
  pipeline = subcluster_pipeline,
  object = sc_object,
  group_col = "cell_type",
  groups = groups_to_run
)

summary_dt <- data.table::rbindlist(purrr::imap(per_group, \(x, name) {
  data.table(
    cell_type = name,
    n_cells = dim(x)[1],
    n_subclusters = length(unique(get_sc_obs(x)$subcluster)),
    hvg_shared_with_global = length(intersect(get_hvg(x), global_hvg))
  )
}))

summary_dt[order(-n_cells)]
#>    cell_type n_cells n_subclusters hvg_shared_with_global
#>       <char>   <int>         <int>                  <int>
#> 1:   T cells    1043             2                   1021
#> 2: Monocytes     459             3                    804
#> 3:        NK     382             3                    697
#> 4:   B cells     279             1                    585
```

Every lineage picked a different feature set, and the shared fraction
with the global HVGs is nowhere near complete. Running one global PCA
and hoping it resolves everything was never going to work.

Writing all of that back is one call.
[`merge_subset_obs()`](https://gregorlueg.github.io/bixverse/reference/merge_subset_obs.md)
takes the list straight from
[`apply_pipeline_per_group()`](https://gregorlueg.github.io/bixverse/reference/apply_pipeline_per_group.md),
checks that no cell shows up twice, and does a single join.

``` r

sc_object <- merge_subset_obs(
  object = sc_object,
  subsets = per_group,
  cols = "subcluster",
  new_names = "fine_label",
  prefix_values = TRUE
)
#> Merged 1 column(s) for 2163 cells into obs (537 cells left as NA).

get_sc_obs(sc_object, filtered = TRUE)[, .N, by = fine_label][order(fine_label)]
#>     fine_label     N
#>         <char> <int>
#> 1:   B cells_0   279
#> 2: Monocytes_0   239
#> 3: Monocytes_1   128
#> 4: Monocytes_2    92
#> 5:        NK_0   140
#> 6:        NK_1   131
#> 7:        NK_2   111
#> 8:   T cells_0   562
#> 9:   T cells_1   481
```

The prefix is doing real work here. Every group’s Leiden run starts
counting at zero, so without it the T cell `1` and the B cell `1`
collapse into one label. Leave it out and you get a warning naming the
values that clash.

``` r

embedding_plot_sc(
  sc_object,
  embedding = "umap",
  colour_by = "fine_label",
  discrete = TRUE
) +
  labs(colour = "Fine label")
```

![](single_cell_subsets_files/figure-html/fine%20umap-1.png)

## Caveats and what is next

A few sharp edges to be aware of before you build on this.

Metacell generation works on subsets too, so
[`generate_bt_meta_cells_sc()`](https://gregorlueg.github.io/bixverse/reference/generate_bt_meta_cells_sc.md),
[`generate_seacells_sc()`](https://gregorlueg.github.io/bixverse/reference/generate_seacells_sc.md)
and
[`generate_supercells_sc()`](https://gregorlueg.github.io/bixverse/reference/generate_supercells_sc.md)
can all be pointed straight at a `SingleCellsSubset`, as can
[`get_pseudobulked_sc()`](https://gregorlueg.github.io/bixverse/reference/get_pseudobulked_sc.md).
Combined with
[`apply_pipeline_per_group()`](https://gregorlueg.github.io/bixverse/reference/apply_pipeline_per_group.md)
that gives you per-cell-type metacells in a few lines. Note that meta
cell memberships come back as positions in the *parent* obs table, which
is what makes them joinable back onto the original object.

Subsets are RNA only. If you hand `modality = "adt"` to a getter you get
an error, so multi-modal work stays on
[`SingleCellsMultiModal`](https://gregorlueg.github.io/bixverse/articles/multi_modal_single_cells.html).

Subsets nest conceptually but not in code: you cannot build a
`SingleCellsSubset` from another `SingleCellsSubset`. If you want a
third level, write the labels back to the parent and subset from there.

## Clean up

``` r

unlink(tempdir_pbmc, recursive = TRUE, force = TRUE)
```
