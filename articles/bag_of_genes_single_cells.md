# Gene set analysis and GRN inference on PBMCs

## Intro

This vignette demonstrates the gene set scoring, spatial
autocorrelation, and gene regulatory network (GRN) inference methods
available in `bixverse`. It picks up where the [PBMC
vignette](https://gregorlueg.github.io/bixverse/articles/single_cell_pbmc.html)
left off: QC, normalisation, PCA, nearest neighbours, clustering and
UMAP are all assumed to be done already. If you have not read the
[design
choices](https://gregorlueg.github.io/bixverse/articles/design_single_cell.html)
and the [introductory
vignette](https://gregorlueg.github.io/bixverse/articles/thinking_single_cell.html),
please do so first.

The methods fall into two broad categories:

- The first, module scores, AUCell, and VISION, takes *pre-defined* gene
  sets and scores them per cell.
- The second, Hotspot and SCENIC, *discovers* gene programmes and
  regulatory relationships from the data itself.
- The third, NMF, also *discovers* gene programmes and regulatory
  relationships via matrix factorisation which is different from SCENIC
  and Hotspot.

We use the PBMC3k dataset throughout. At 2,700 cells this is likely too
small for the GRN methods to produce biologically meaningful results,
but it is large enough to show the API and the workflow end to end. On
real datasets with tens of thousands of cells the same code applies
unchanged. Also, as the

``` r

library(bixverse)
library(ggplot2)
library(data.table)
#> 
#> Attaching package: 'data.table'
#> The following object is masked from 'package:base':
#> 
#>     %notin%
library(magrittr)
```

> **Note**
>
> The vignette was built on a GitHub runner which is a 2-core machine
> with ~8 GB of memory. That this works is incredible enough, but take
> this into account when looking at this vignette.

### Rebuilding the processed object

We reconstruct the processed `SingleCells` object from the PBMC
vignette. The chunk is hidden for brevity: it downloads the data, runs
QC, HVG selection, PCA, neighbour computation, Leiden clustering and
UMAP.

Rebuild the PBMC3k object (click to expand)

``` r

pbmc3k_path <- download_pbmc3k()
tempdir_pbmc <- tempdir()

sc_object <- SingleCells(dir_data = tempdir_pbmc)
mtx_io_params <- get_cell_ranger_params(pbmc3k_path)

sc_object <- load_mtx(
  object = sc_object,
  sc_mtx_io_param = mtx_io_params,
  mtx_streaming = FALSE,
  .verbose = FALSE
)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

setnames_sc(sc_object, table = "var", old = "column1", new = "gene_symbol")
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

var <- get_sc_var(sc_object)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
ensembl_to_symbol <- setNames(var$gene_symbol, var$gene_id)
symbol_to_ensembl <- setNames(var$gene_id, var$gene_symbol)

# gene set proportions for QC
gs_of_interest <- list(
  MT = var[grepl("^MT-", gene_symbol), gene_id],
  Ribo = var[grepl("^RPS|^RPL", gene_symbol), gene_id]
)
sc_object <- gene_set_proportions_sc(
  sc_object,
  gs_of_interest,
  streaming = FALSE,
  .verbose = FALSE
)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

# MAD QC
qc_df <- sc_object[[c("cell_id", "lib_size", "nnz", "MT")]]
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
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
qc <- run_cell_qc(
  metrics = metrics,
  cells_to_keep = get_cells_to_keep(sc_object),
  directions = directions,
  threshold = 3
)
sc_object[["outlier"]] <- qc$combined
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
cells_to_keep <- qc_df[!qc$combined, cell_id]
sc_object <- set_cells_to_keep(sc_object, cells_to_keep)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

# HVG, PCA, neighbours, clustering, UMAP
sc_object <- find_hvg_sc(sc_object, hvg_no = 2000L, .verbose = FALSE)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
sc_object <- calculate_pca_sc(sc_object, no_pcs = 30L, sparse_svd = TRUE)
#> Using sparse SVD solving on scaled data on 2000 HVG.
sc_object <- find_neighbours_sc(
  sc_object,
  neighbours_params = params_sc_neighbours(
    knn = list(knn_method = "exhaustive")
  )
)
#> 
#> Generating sNN graph (full: TRUE).
#> Transforming sNN data to igraph.
sc_object <- find_clusters_sc(sc_object, res = 1.5, name = "leiden_clusters")
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
sc_object <- umap_sc(sc_object, knn_method = "annoy")
#> Running UMAP.
#> Using n_epochs = 500 (dataset <10k samples or adam_parallel optimiser)
#> Using provided kNN graph.

# UMAP coordinates for reuse
umap_dt <- extract_embedding_data(sc_object, embedding = "umap")
setnames(
  umap_dt,
  old = c("dim_1", "dim_2"),
  new = c("umap_1", "umap_2"),
  skip_absent = TRUE
)
```

## Gene set scoring

### Module scores

The simplest approach to gene set activity is the module score from
[Tirosh et
al. (2016)](https://www.science.org/doi/10.1126/science.aad0501): for
each cell, compute the mean expression of the gene set and subtract the
mean expression of a size-matched control set drawn from the same
expression bins. We will just use some lineage-specific genes to have
pretty visualisations.

``` r

lineage_sets <- list(
  `T cell` = symbol_to_ensembl[c(
    "CD3D",
    "CD3E",
    "CD3G",
    "CD2",
    "IL7R",
    "CD7",
    "LEF1",
    "TCF7",
    "LTB",
    "TRAC"
  )],
  `Cytotoxic NK` = symbol_to_ensembl[c(
    "NKG7",
    "GNLY",
    "GZMA",
    "GZMB",
    "GZMH",
    "PRF1",
    "CST7",
    "KLRD1",
    "KLRB1",
    "FGFBP2"
  )],
  `B cell` = symbol_to_ensembl[c(
    "CD79A",
    "CD79B",
    "MS4A1",
    "BANK1",
    "IGHM",
    "IGHD",
    "CD74",
    "HLA-DQA1",
    "TCL1A",
    "VPREB3"
  )],
  `Monocyte` = symbol_to_ensembl[c(
    "CD14",
    "LYZ",
    "S100A8",
    "S100A9",
    "S100A12",
    "CST3",
    "FCN1",
    "VCAN",
    "MNDA",
    "TYROBP"
  )]
)

# drop any NAs
lineage_sets <- lapply(lineage_sets, function(x) x[!is.na(x)])

module_scores <- module_scores_sc(
  object = sc_object,
  gs_list = lineage_sets,
  .verbose = FALSE
)

# you need to unclass here
ms_dt <- as.data.table(unclass(module_scores), keep.rownames = "cell_id")
ms_dt <- merge(ms_dt, umap_dt, by = "cell_id")

ms_long <- melt(
  ms_dt,
  id.vars = c("cell_id", "umap_1", "umap_2"),
  measure.vars = names(lineage_sets),
  variable.name = "phase",
  value.name = "score"
)

# scale the scores to make them more comparable across
ms_long[,
  score_scaled := (score - min(score)) / (max(score) - min(score)),
  by = phase
]

ggplot(ms_long, aes(x = umap_1, y = umap_2)) +
  geom_point(aes(colour = score_scaled), size = 0.3) +
  scale_colour_viridis_c() +
  facet_wrap(~phase, ncol = 2) +
  theme_bw() +
  labs(colour = "Module\nscore\n(scaled)")
```

![](bag_of_genes_single_cells_files/figure-html/module-scores-1.png)

Module scores for lineage genes projected onto UMAP

Let’s imagine you want to add the scores to the internal DB. There are
helpers for this, too. The
[`get_data()`](https://gregorlueg.github.io/bixverse/reference/get_data.md)
helper in combination with \`\`

``` r

sc_object <- add_sc_new_obs(
  object = sc_object,
  obs_data = get_data(module_scores)
)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

head(sc_object)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#>    cell_idx          cell_id   nnz lib_size to_keep          MT      Ribo
#>       <num>           <char> <num>    <num>  <lgcl>       <num>     <num>
#> 1:        1 AAACATACAACCAC-1   778     2418    TRUE 0.030190241 0.4371381
#> 2:        3 AAACATTGATCAGC-1  1126     3144    TRUE 0.008905852 0.3171120
#> 3:        4 AAACCGTGCTTCCG-1   953     2632    TRUE 0.017477203 0.2431611
#> 4:        6 AAACGCACTGGTAC-1   779     2154    TRUE 0.016713092 0.3635097
#> 5:        8 AAACGCTGGTTCTT-1   785     2255    TRUE 0.031042129 0.3844789
#> 6:        9 AAACGCTGTAGCCA-1   530     1273    TRUE 0.011783189 0.3794187
#>    outlier leiden_clusters     T cell Cytotoxic NK     B cell   Monocyte
#>     <lgcl>           <int>      <num>        <num>      <num>      <num>
#> 1:   FALSE               1  0.3236977   -0.2325926 -0.3144407 -0.7570705
#> 2:   FALSE               1  0.2307960   -0.3455873 -0.4134143 -0.6403040
#> 3:   FALSE               7 -1.1175209   -0.2854205  0.0692758  1.2887695
#> 4:   FALSE               1 -0.1862041   -0.5217526 -0.2286735 -0.7023786
#> 5:   FALSE               4  0.1805688    0.4830460 -0.4746688 -0.6906318
#> 6:   FALSE               4 -0.1726186    0.2581931  0.7847812 -0.7003754
```

As we can appreciate, this has now added the new columns.

### AUCell

[AUCell](https://www.nature.com/articles/nmeth.4463) ranks genes within
each cell by expression and summarises where the gene set lands in that
ranking. More robust to outliers than a simple mean. `bixverse` gives
you three statistics via
[`params_sc_aucell()`](https://gregorlueg.github.io/bixverse/reference/params_sc_aucell.md),
all reading the same ranking but weighting it differently:

- `"recovery"`, the default, is the recovery-curve AUC under a rank
  cutoff, i.e. what Aibar et al. published. Top-heavy: only genes inside
  the top `max_rank` contribute, so it is the sharper instrument for
  marker sets and the one to use for SCENIC.
- `"wilcox"` is the AUC from the Mann-Whitney U statistic over the full
  ranking. Its null sits at 0.5 whatever the gene set size, which makes
  it a safe pick for pathway activity, but the flatter score does not
  split into on/off populations.
- `"ap"` is average precision. Most top-heavy of the three, but its null
  tracks gene set prevalence, so don’t compare raw values across sets of
  different size unless you switch `standardise` on.

We stick with the default here.

``` r

auc_scores <- aucell_sc(
  object = sc_object,
  gs_list = lineage_sets,
  aucell_params = params_sc_aucell(),
  .verbose = FALSE
)

auc_dt <- as.data.table(unclass(auc_scores), keep.rownames = "cell_id")
auc_dt <- merge(auc_dt, umap_dt, by = "cell_id")

auc_long <- melt(
  auc_dt,
  id.vars = c("cell_id", "umap_1", "umap_2"),
  measure.vars = names(lineage_sets),
  variable.name = "phase",
  value.name = "auc"
)

ggplot(auc_long, aes(x = umap_1, y = umap_2)) +
  geom_point(aes(colour = auc), size = 0.3) +
  scale_colour_viridis_c() +
  facet_wrap(~phase) +
  theme_bw() +
  labs(colour = "AUC")
```

![](bag_of_genes_single_cells_files/figure-html/aucell-1.png)

AUCell scores for lineage genes projected onto UMAP

### VISION with spatial autocorrelation

Scoring gene sets is useful, but a natural follow-up question is: *are
these scores spatially structured on the cell neighbourhood graph, or
just noise?*
[VISION](https://www.nature.com/articles/s41467-019-12235-0) provides a
permutation-based test of spatial autocorrelation (Geary’s C) on the kNN
graph. Gene sets with significant autocorrelation are those whose
activity varies smoothly across cell neighbourhoods – a strong
indication that the signal is biologically meaningful rather than
random.

VISION expects gene sets in a signed format with `pos` and (optionally)
`neg` components. In the cause of our cell lineage genes, we could use
the other cell type’s markers as negatives, but we will just leave it
blank here.

``` r

vision_gs <- lapply(lineage_sets, function(genes) list(pos = genes))

vision_res <- vision_w_autocor_sc(
  object = sc_object,
  gs_list = vision_gs,
  embd_to_use = "pca",
  vision_params = params_sc_vision(n_perm = 500L),
  .verbose = FALSE
)

head(vision_res$auto_cor_dt)
#>    gene_set_name  auto_cor       p_val         fdr
#>           <char>     <num>       <num>       <num>
#> 1:        T cell 0.7612149 0.001996008 0.003992016
#> 2:  Cytotoxic NK 0.4620488 0.009980040 0.009980040
#> 3:        B cell 0.5621615 0.009980040 0.009980040
#> 4:      Monocyte 0.5473183 0.001996008 0.003992016
```

Unsurprisingly, all of these gene sets show highly significant spatial
correlation.

## Hotspot

Where the methods above test *pre-defined* gene sets, Hotspot ([DeTomaso
& Yosef,
2021](https://www.cell.com/cell-systems/fulltext/S2405-4712(21)00114-9))
discovers gene programmes *de novo* by testing each gene for spatial
autocorrelation on the kNN graph and then grouping the significant genes
by their local correlation structure.

### Gene autocorrelation

The first step computes Geary’s C for every gene against the neighbour
graph. Genes with significant autocorrelation are those whose expression
varies smoothly across neighbouring cells, i.e., genes that mark
spatially coherent programmes.

``` r

hotspot_autocor <- hotspot_autocor_sc(
  object = sc_object,
  .verbose = FALSE
)

hotspot_autocor[, gene_symbol := ensembl_to_symbol[gene_id]]

head(hotspot_autocor[order(fdr)], 25L)
#>             gene_id  gaerys_c   z_score  pval   fdr gene_symbol
#>              <char>     <num>     <num> <num> <num>      <char>
#>  1: ENSG00000188290 0.3242266  39.66046     0     0        HES4
#>  2: ENSG00000119535 0.3022016  49.27494     0     0       CSF3R
#>  3: ENSG00000163131 0.5119517  63.60383     0     0        CTSS
#>  4: ENSG00000163191 0.4122631  58.87283     0     0     S100A11
#>  5: ENSG00000163220 0.6936087 109.61599     0     0      S100A9
#>  6: ENSG00000143546 0.6563530 114.56873     0     0      S100A8
#>  7: ENSG00000197956 0.6433797  98.70326     0     0      S100A6
#>  8: ENSG00000196154 0.6429363  92.52890     0     0      S100A4
#>  9: ENSG00000177954 0.6083654  94.29710     0     0       RPS27
#> 10: ENSG00000158869 0.6997935  81.83833     0     0      FCER1G
#> 11: ENSG00000203747 0.6159353  71.61398     0     0      FCGR3A
#> 12: ENSG00000198821 0.2228819  59.91400     0     0       CD247
#> 13: ENSG00000143185 0.3826069  53.56304     0     0        XCL2
#> 14: ENSG00000143184 0.3719576  57.30895     0     0        XCL1
#> 15: ENSG00000116667 0.2143988  43.89426     0     0     C1orf21
#> 16: ENSG00000143947 0.3742483  57.71536     0     0      RPS27A
#> 17: ENSG00000115523 0.5908421 126.38232     0     0        GNLY
#> 18: ENSG00000153563 0.2527090  47.25644     0     0        CD8A
#> 19: ENSG00000172116 0.2330091  39.71748     0     0        CD8B
#> 20: ENSG00000071082 0.4138269  64.81413     0     0       RPL31
#> 21: ENSG00000144713 0.3622806  61.70233     0     0       RPL32
#> 22: ENSG00000168028 0.2855044  45.84666     0     0        RPSA
#> 23: ENSG00000233276 0.4964235  64.97836     0     0        GPX1
#> 24: ENSG00000163931 0.2321878  38.70044     0     0         TKT
#> 25: ENSG00000196542 0.3874654  51.78241     0     0      SPTSSB
#>             gene_id  gaerys_c   z_score  pval   fdr gene_symbol
#>              <char>     <num>     <num> <num> <num>      <char>
```

You should get a set of well known cell type markers here…

### Gene-gene local correlations and modules

Taking the significant genes forward, we compute pairwise local
correlations and cluster them into modules. Each module represents a
group of genes that are not only individually autocorrelated but also
co-vary locally - a much stronger signal than global correlation alone.

``` r

sig_genes <- hotspot_autocor[fdr <= 0.05, gene_id]

hotspot_cor <- hotspot_gene_cor_sc(
  object = sc_object,
  genes_to_take = sig_genes,
  .verbose = TRUE
)

hotspot_cor
#> Hotspot gene-gene local correlation results
#>   Genes: 2320
#>   Cells: 2163
#>   Modules: not yet computed (see generate_hotspot_membership)
```

This returns a HotSpot result. Let’s add the membership and plot it.

``` r

hotspot_cor <- generate_hotspot_membership(hotspot_cor)

# this will only plot a subsample of 500 genes for speed (stratified by
# module). You can control this via a parameter in the plotting function
plot(hotspot_cor)
```

![](bag_of_genes_single_cells_files/figure-html/hotspot%20-%20plot%20modules-1.png)

Let’s pull out the genes:

``` r

membership <- get_hotspot_membership(hotspot_cor)
membership[, gene_symbol := ensembl_to_symbol[gene_id]]

membership <- membership[!is.na(cluster_member)]

head(membership[order(cluster_member)], 20L)
#>             gene_id cluster_member gene_symbol
#>              <char>          <num>      <char>
#>  1: ENSG00000137959            477      IFI44L
#>  2: ENSG00000134321            477       RSAD2
#>  3: ENSG00000115415            477       STAT1
#>  4: ENSG00000196141            477     SPATS2L
#>  5: ENSG00000136514            477        RTP4
#>  6: ENSG00000138642            477       HERC6
#>  7: ENSG00000138646            477       HERC5
#>  8: ENSG00000003147            477        ICA1
#>  9: ENSG00000119917            477       IFIT3
#> 10: ENSG00000185745            477       IFIT1
#> 11: ENSG00000173456            477       RNF26
#> 12: ENSG00000135114            477        OASL
#> 13: ENSG00000140511            477      HAPLN3
#> 14: ENSG00000157601            477         MX1
#> 15: ENSG00000181817            515       LSM10
#> 16: ENSG00000173193            515      PARP14
#> 17: ENSG00000113719            515      ERGIC1
#> 18: ENSG00000213722            515       DDAH2
#> 19: ENSG00000146192            515        FGD2
#> 20: ENSG00000112335            515        SNX3
#>             gene_id cluster_member gene_symbol
#>              <char>          <num>      <char>
```

The resulting modules should roughly correspond to the major cell-type
programmes in PBMCs (T cell, monocyte, B cell, NK cell signatures,
etc.). Let’s visualise this with AUCell check it out.

``` r

hotspot_gene_sets <- membership %$% split(gene_id, cluster_member)

hotspot_gene_sets <- lapply(hotspot_gene_sets, function(genes) {
  list(pos = genes)
})

vision_scores_hotspot <- vision_sc(
  object = sc_object,
  gs_list = hotspot_gene_sets,
  .verbose = FALSE
)

vision_scores_dt <- as.data.table(
  unclass(vision_scores_hotspot),
  keep.rownames = "cell_id"
)
vision_scores_dt <- merge(vision_scores_dt, umap_dt, by = "cell_id")


vision_scores_long <- melt(
  vision_scores_dt,
  id.vars = c("cell_id", "umap_1", "umap_2"),
  measure.vars = names(hotspot_gene_sets),
  variable.name = "gene_set",
  value.name = "vision_score"
)

# min max individually for prettier visualisations
vision_scores_long[,
  vision_score_scaled := (vision_score - min(vision_score)) /
    (max(vision_score) - min(vision_score)),
  by = gene_set
]

ggplot(vision_scores_long, aes(x = umap_1, y = umap_2)) +
  geom_point(aes(colour = vision_score_scaled), size = 0.3) +
  scale_colour_viridis_c() +
  facet_wrap(~gene_set) +
  theme_bw() +
  labs(colour = "Vision (scaled)")
```

![](bag_of_genes_single_cells_files/figure-html/hotspot%20-%20plot%20the%20gene%20sets-1.png)

## SCENIC

[SCENIC](https://www.nature.com/articles/nmeth.4463) infers gene
regulatory networks by asking: for each gene, which transcription
factors (TFs) best predict its expression? The original implementation
uses either random forests or GRNBoost2 (gradient boosted trees) to
regress each target gene on all TFs, then extracts feature importances
as a proxy for regulatory strength.

The `bixverse` implementation re-implements this in Rust with several
optimisations that make the RF and ExtraTrees paths substantially faster
than the single-target-at-a-time approach used by GRNBoost2:

- **Quantisation**: TF expression values are discretised into 256 bins
  (u8), so histogram construction during tree building operates on
  single bytes and fits comfortably in L1 cache.
- **Multi-output batching**: up to 64 target genes share a single tree
  structure. The histogram construction cost is paid once per node per
  feature, but the split score aggregates variance reduction across all
  targets in the batch. This amortises the dominant cost of tree
  building across targets.
- **Correlated gene batching**: rather than assigning genes to batches
  randomly, an optional SVD + k-means step groups co-expressed genes
  together so that the shared tree structure is more informative for
  each target in the batch.
- **GRNBoost2 with histogram subtraction**: for the gradient boosted
  path (single-target), the code builds full-feature histograms at each
  node and derives the larger child via parent-minus-smaller
  subtraction, avoiding a redundant scan of the larger child’s samples.
  OOB early stopping prevents overfitting and is the main source of
  speedup.

On the 2,700-cell PBMC3k dataset these optimisations are not that
noticeable… The data is too small for any method to take long. On
datasets with tens of thousands of cells and thousands of TFs, the RF/ET
multi-output path comfortably outperforms GRNBoost2. You need to just
decide if you are okay with the batching of genes. Generally speaking,
the big signals will be recovered again and again (from empirical
testing).

### Gene filtering

SCENIC first filters genes by minimum total counts and minimum
proportion of expressing cells to remove uninformative targets.

``` r

scenic_genes <- scenic_gene_filter_sc(
  object = sc_object,
  scenic_params = params_scenic(min_counts = 50L),
  .verbose = FALSE
)

length(scenic_genes)
#> [1] 5430
```

### Transcription factor list

We need a list of known TFs. The Aerts lab provides a curated list for
human which we can download and map to Ensembl IDs.

``` r

tf_dt <- data.table::fread(
  "https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt",
  header = FALSE,
  col.names = "tf"
)
tf_dt[, gene_id := symbol_to_ensembl[tf]]
tf_dt <- tf_dt[!is.na(gene_id)]

nrow(tf_dt)
#> [1] 1100
```

### GRN inference

With genes and TFs in hand, we run the GRN inference step. Here we use
random forests with a batch size of 64 to illustrate the multi-output
path.

``` r

scenic_res <- scenic_grn_sc(
  object = sc_object,
  tf_ids = tf_dt$gene_id,
  genes_to_take = scenic_genes,
  scenic_params = params_scenic(
    learner_type = "randomforest",
    gene_batch_size = 64L
  ),
  .verbose = TRUE
)
#> SCENIC: 5430 target genes, 466 TFs, 2163 cells

scenic_res
#> ScenicGrn (GRN results)
#>   No genes:                 5430 
#>   No TFs:                   466 
#>   Applied learner:          randomforest 
#>   TF to gene generated:     FALSE 
#>   CisTarget res generated:  FALSE
```

If you are bored, you can run both methods
(`learner_type = "grnboost2"`) and compare the importance scores. You
will see very high correlations here (and can explore the speed
differences…)

### TF-to-gene refinement

The raw importance matrix is large and noisy. The next steps winnow it
down to a candidate regulatory network by keeping only the top TF-gene
pairs and filtering by pairwise expression correlation. The approach
used here is that per gene only TFs with an
`importance score ≥ mean(importance_score) + n_sd * sd(importance_score)`
per given gene are kept. Note that `n_sd` should be adjusted depending
on the learner: RF and ExtraTrees spread importance more evenly across
TFs than GRNBoost2, so a lower threshold (e.g. `n_sd = 1`) is
appropriate, whereas GRNBoost2 concentrates importance onto fewer TFs
and tolerates `n_sd = 1.5` or even `n_sd = 2`.

``` r

scenic_res <- identify_tf_to_genes(
  scenic_res,
  n_sd = 2,
  .verbose = TRUE
)
#> Extracting TF to gene associations via per-gene threshold (mean + 2.0 * SD).

scenic_res <- tf_to_genes_correlations(
  x = scenic_res,
  object = sc_object,
  .verbose = TRUE
)
#> Calculating the pairwise correlations between the TFs and genes
#> Keeping activating TF <> gene links at |rho| > 0.030
#> Removing self loops (TF controlling its own expression

tf_gene_dt <- get_tf_to_gene(scenic_res)
tf_gene_dt[, tf_symbol := ensembl_to_symbol[tf]]
tf_gene_dt[, gene_symbol := ensembl_to_symbol[gene]]

head(tf_gene_dt[order(-importance)], 5L)
#>                 tf            gene importance pairwise_cor cor_sign tf_symbol
#>             <char>          <char>      <num>        <num>    <int>    <char>
#> 1: ENSG00000171223 ENSG00000120129  0.2948482    0.3103157        1      JUNB
#> 2: ENSG00000066336 ENSG00000107341  0.2746513    0.1931868        1      SPI1
#> 3: ENSG00000170345 ENSG00000120129  0.2573821    0.3420300        1       FOS
#> 4: ENSG00000066336 ENSG00000165025  0.2460342    0.2060715        1      SPI1
#> 5: ENSG00000139187 ENSG00000153563  0.2453607    0.2005602        1     KLRG1
#>    gene_symbol
#>         <char>
#> 1:       DUSP1
#> 2:      UBE2R2
#> 3:       DUSP1
#> 4:         SYK
#> 5:        CD8A
```

### CisTarget motif enrichment

The final (and optional) step in the SCENIC pipeline validates the
predicted TF-gene links against known motif binding sites. This requires
downloading the CisTarget motif ranking database and motif-to-TF
annotation files from the [Aerts lab
resources](https://resources.aertslab.org/cistarget/), which together
are roughly 400 MB. We therefore show the code but do not evaluate it by
default.

``` r

# Download the CisTarget reference files (cached after first download)
paths <- download_cistarget_hg38()

rankings <- read_motif_ranking(paths$rankings)
annotations <- read_motif_annotation_file(paths$motif_annotations)

# Run motif enrichment on the TF-gene links
scenic_res <- tf_to_genes_motif_enrichment(
  x = scenic_res,
  motif_rankings = rankings,
  annot_data = annotations,
  # we store everything as Ensembl... The SCENIC data has gene symbols, hence,
  # the need to provide conversion
  gene_id_to_symbol = ensembl_to_symbol,
  cis_target_params = params_cistarget(),
  only_high_conf_tf = FALSE,
  .verbose = TRUE
)

# CisTarget results per TF
cistarget_dt <- get_cistarget_res(scenic_res)
head(cistarget_dt[order(-nes)])

# The TF-to-gene table now has a leading_edge column indicating which
# genes fall within the motif recovery curve
tf_gene_final <- get_tf_to_gene(scenic_res)

head(tf_gene_final[(in_leading_edge)][order(-importance)])
```

CisTarget is not a pass/fail filter on the module. Per enriched motif it
builds a recovery curve, compares it against the mean + 2 SD curve
across every motif in the database, and keeps only the genes above the
point of maximum separation. So expect the regulons coming out of
[`build_regulons()`](https://gregorlueg.github.io/bixverse/reference/build_regulons.md)
to be a good deal smaller than the modules going in. This is how you
would run the regulon binarisation.

``` r

# build_regulons() applies that filter, adds the TF back to its own set and
# drops anything too small. This is what you hand to aucell_sc()
regulons <- build_regulons(scenic_res)

auc_regulons <- aucell_sc(
  object = sc_object,
  gs_list = regulons,
  aucell_params = params_sc_aucell()
)

binarised_regulons <- binarise_regulon_activity(auc_regulons)

binary_matrix <- binarised_regulons$binary
colnames(binary_matrix) <- ensembl_to_symbol[colnames(binary_matrix)]

sc_object[["STAT1_active"]] <- binary_matrix[, "STAT1"]

embedding_plot_sc(
  object = sc_object,
  embedding = "umap",
  colour_by = "STAT1_active"
)
```

This is how you work with all types of `"bag of genes"` analyses for
single cell in `bixverse`.

## NMF on single cells

The methods so far have all started from either a pre-defined gene set,
a graph-based notion of locality, or a TF-centric regression framework.
NMF (non-negative matrix factorisation) takes yet another angle: it
factorises the expression matrix `V ≈ W H` into a small number of latent
factors, where `W` carries non-negative gene loadings and `H` carries
non-negative cell activations. The non-negativity constraint tends to
push the factors into something interpretable as additive gene
programmes rather than abstract directions of variance like PCA.
`bixverse` implements
[HALS](https://proceedings.mlr.press/v39/kimura14.pdf) (Hierarchical
Alternating Least Squares) in Rust, which converges fast and scales
nicely.

A few practical notes before diving in:

- NMF is non-convex: different random initialisations land in different
  local optima. The package therefore provides
  [`nmf_sc()`](https://gregorlueg.github.io/bixverse/reference/nmf_sc.md)
  for a single run and
  [`stabilised_nmf_sc()`](https://gregorlueg.github.io/bixverse/reference/stabilised_nmf_sc.md)
  for multiple restarts. The latter is the right default for any serious
  analysis as well as down stream meta analyses such as finding
  consensus programmes.
- Picking `k` is the eternal question. Too small and programmes get
  smushed together; too large and they fragment into redundant copies.
  There is no free lunch here. Usually you run with several values and
  inspect. If in doubt, lick your finger and see from where the wind
  comes… (Joke)
- NMF on a heterogeneous dataset often just rediscovers the obvious
  cell-type axis. The interesting use case is running it *within* a cell
  type to find finer-grained programmes (activation, exhaustion,
  cycling, etc.). We will do exactly that here.

We will focus on the T cells as an example. Think of your own use cases
in your data set and analyses.

### Identifying T cells

We will use the AUC scores from earlier to pull out the T cell
population. [Otsu’s
method](https://en.wikipedia.org/wiki/Otsu%27s_method) gives us a
data-driven threshold for the bimodal AUC distribution.

``` r

otsu_t_cell <- find_threshold_otsu(auc_dt$`T cell`)

ggplot(data = auc_dt, mapping = aes(x = `T cell`)) +
  geom_histogram(bins = 13, fill = "darkgrey") +
  theme_minimal() +
  xlab("T cell signature [AUC]") +
  geom_vline(xintercept = otsu_t_cell)
```

![](bag_of_genes_single_cells_files/figure-html/define%20t-cell%20threshold-1.png)

With the threshold in hand we grab the T cell IDs and recompute HVGs
*within* this subset. This matters: HVGs computed across all PBMCs are
dominated by inter-cell-type variation (T vs B vs monocyte). HVGs within
T cells expose the genes that actually vary across T cell states, which
is what we want NMF to factorise.

[`get_hvg_data_sc()`](https://gregorlueg.github.io/bixverse/reference/get_hvg_data_sc.md)
is a non-mutating variant of
[`find_hvg_sc()`](https://gregorlueg.github.io/bixverse/reference/find_hvg_sc.md):
it returns a data.table with the HVG metrics and an `is_hvg` flag,
without overwriting the HVGs stored on the object. Convenient when you
want HVGs for a specific downstream analysis without disturbing the
global state.

``` r

t_cell_ids <- auc_dt[`T cell` >= otsu_t_cell, cell_id]

hvg_t_cell_data <- get_hvg_data_sc(object = sc_object, cell_ids = t_cell_ids)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpCS7SaE/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

t_cell_hvg <- hvg_t_cell_data[(is_hvg), gene_id]
```

### Running a single NMF run

[`nmf_sc()`](https://gregorlueg.github.io/bixverse/reference/nmf_sc.md)
runs one NMF on the chosen cells and genes. The default initialisation
is NNDSVD, which is deterministic and usually gives a sensible starting
point.

``` r

t_cell_nmf_results <- nmf_sc(
  object = sc_object,
  k = 10L,
  cell_ids = t_cell_ids,
  gene_ids = t_cell_hvg
)
```

The result is an `NmfResult` carrying the `W` gene-loading matrix (genes
x components) and the `H` cell-activation matrix (components x cells),
plus the parameters and convergence info.

``` r

get_w(t_cell_nmf_results)[1:5, 1:5]
#>                      comp_01      comp_02   comp_03      comp_04      comp_05
#> ENSG00000188976 2.558843e+00 9.992267e-11 0.3763380 2.132400e-01 9.993509e-11
#> ENSG00000187608 9.983787e-11 9.992267e-11 0.5599045 9.981262e-11 9.993509e-11
#> ENSG00000186827 9.983787e-11 9.992267e-11 7.3252878 9.981262e-11 9.993509e-11
#> ENSG00000176022 1.122509e-01 2.855008e-02 0.1625732 9.981262e-11 9.993509e-11
#> ENSG00000242485 2.209524e+00 2.744842e+00 4.3223567 6.435678e-01 1.557151e+00
```

``` r

get_h(t_cell_nmf_results)[1:5, 1:5]
#>         AAACATACAACCAC-1 AAACATTGATCAGC-1 AAACGCTGGTTCTT-1 AAACTTGATCCAGA-1
#> comp_01       0.03039108      0.028655788     2.779504e-02     0.0335300341
#> comp_02       0.04223846      0.009922862     6.078151e-02     0.0007059687
#> comp_03       0.01914044      0.041368406     1.002924e-10     0.0018196171
#> comp_04       0.04766104      0.033128884     3.672349e-02     0.0398822948
#> comp_05       0.03773499      0.001358267     4.147765e-02     0.0049292077
#>         AAAGCCTGTATGCG-1
#> comp_01       0.03014318
#> comp_02       0.02215907
#> comp_03       0.04835882
#> comp_04       0.02421124
#> comp_05       0.04118023
```

### Running multiple NMF runs

For anything beyond exploratory work you want multiple restarts.
[`stabilised_nmf_sc()`](https://gregorlueg.github.io/bixverse/reference/stabilised_nmf_sc.md)
runs `n_runs` HALS NMFs with random initialisations (seeded by
`seed + i`) and column-binds the resulting `W` matrices. Looking at this
stabilised data allows you to identify components that constantly show
up between random initialisations.

``` r

t_cell_nmf_results_stabilised <- stabilised_nmf_sc(
  object = sc_object,
  k = 10L,
  cell_ids = t_cell_ids,
  gene_ids = t_cell_hvg,
  n_runs = 10L
)
```

The W matrix here is wide: `n_genes x (k * n_runs)`. Each block of `k`
columns corresponds to one run. You can then use it for downstream
consensus analyses.

``` r

get_w(t_cell_nmf_results_stabilised)[1:5, 1:5]
#>                 run_01.comp_01 run_01.comp_02 run_01.comp_03 run_01.comp_04
#> ENSG00000188976   8.834997e-01   1.168006e+00   4.921125e-01   9.998284e-11
#> ENSG00000187608   9.988405e-11   1.000315e-10   4.467033e-02   9.998284e-11
#> ENSG00000186827   9.988405e-11   1.000315e-10   9.998983e-11   9.998284e-11
#> ENSG00000176022   4.819395e-01   1.243956e-01   9.209826e-02   9.998284e-11
#> ENSG00000242485   7.240750e-01   1.973956e+00   1.215576e+00   6.745355e-01
#>                 run_01.comp_05
#> ENSG00000188976   9.488367e-01
#> ENSG00000187608   1.238252e+00
#> ENSG00000186827   1.001258e-10
#> ENSG00000176022   1.001258e-10
#> ENSG00000242485   1.621210e+00
```

The H matrices are kept per-run, since each random initialisation gives
a different cell-activation pattern (components are not aligned across
runs without further work).

``` r

str(get_h(t_cell_nmf_results_stabilised))
#> List of 10
#>  $ run_01: num [1:10, 1:1034] 0.03552 0.00777 0.01237 0.03586 0.03921 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:10] "comp_01" "comp_02" "comp_03" "comp_04" ...
#>   .. ..$ : chr [1:1034] "AAACATACAACCAC-1" "AAACATTGATCAGC-1" "AAACGCTGGTTCTT-1" "AAACTTGATCCAGA-1" ...
#>  $ run_02: num [1:10, 1:1034] 4.08e-02 4.75e-03 3.12e-02 4.11e-03 1.00e-10 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:10] "comp_01" "comp_02" "comp_03" "comp_04" ...
#>   .. ..$ : chr [1:1034] "AAACATACAACCAC-1" "AAACATTGATCAGC-1" "AAACGCTGGTTCTT-1" "AAACTTGATCCAGA-1" ...
#>  $ run_03: num [1:10, 1:1034] 4.34e-02 2.35e-03 1.00e-10 4.04e-02 1.00e-10 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:10] "comp_01" "comp_02" "comp_03" "comp_04" ...
#>   .. ..$ : chr [1:1034] "AAACATACAACCAC-1" "AAACATTGATCAGC-1" "AAACGCTGGTTCTT-1" "AAACTTGATCCAGA-1" ...
#>  $ run_04: num [1:10, 1:1034] 3.93e-02 1.00e-10 5.72e-02 1.35e-02 1.00e-10 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:10] "comp_01" "comp_02" "comp_03" "comp_04" ...
#>   .. ..$ : chr [1:1034] "AAACATACAACCAC-1" "AAACATTGATCAGC-1" "AAACGCTGGTTCTT-1" "AAACTTGATCCAGA-1" ...
#>  $ run_05: num [1:10, 1:1034] 4.20e-02 9.28e-04 1.00e-10 9.99e-11 2.34e-02 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:10] "comp_01" "comp_02" "comp_03" "comp_04" ...
#>   .. ..$ : chr [1:1034] "AAACATACAACCAC-1" "AAACATTGATCAGC-1" "AAACGCTGGTTCTT-1" "AAACTTGATCCAGA-1" ...
#>  $ run_06: num [1:10, 1:1034] 4.65e-02 1.00e-10 1.00e-10 4.25e-02 1.00e-10 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:10] "comp_01" "comp_02" "comp_03" "comp_04" ...
#>   .. ..$ : chr [1:1034] "AAACATACAACCAC-1" "AAACATTGATCAGC-1" "AAACGCTGGTTCTT-1" "AAACTTGATCCAGA-1" ...
#>  $ run_07: num [1:10, 1:1034] 4.60e-02 1.00e-10 1.00e-10 7.62e-02 4.03e-04 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:10] "comp_01" "comp_02" "comp_03" "comp_04" ...
#>   .. ..$ : chr [1:1034] "AAACATACAACCAC-1" "AAACATTGATCAGC-1" "AAACGCTGGTTCTT-1" "AAACTTGATCCAGA-1" ...
#>  $ run_08: num [1:10, 1:1034] 4.49e-02 2.58e-03 9.99e-11 3.42e-02 1.52e-02 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:10] "comp_01" "comp_02" "comp_03" "comp_04" ...
#>   .. ..$ : chr [1:1034] "AAACATACAACCAC-1" "AAACATTGATCAGC-1" "AAACGCTGGTTCTT-1" "AAACTTGATCCAGA-1" ...
#>  $ run_09: num [1:10, 1:1034] 2.87e-02 2.09e-02 3.28e-02 2.82e-02 1.00e-10 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:10] "comp_01" "comp_02" "comp_03" "comp_04" ...
#>   .. ..$ : chr [1:1034] "AAACATACAACCAC-1" "AAACATTGATCAGC-1" "AAACGCTGGTTCTT-1" "AAACTTGATCCAGA-1" ...
#>  $ run_10: num [1:10, 1:1034] 0.0102 0.0463 0.0281 0.0379 0.0254 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:10] "comp_01" "comp_02" "comp_03" "comp_04" ...
#>   .. ..$ : chr [1:1034] "AAACATACAACCAC-1" "AAACATTGATCAGC-1" "AAACGCTGGTTCTT-1" "AAACTTGATCCAGA-1" ...
```

The run with the lowest final reconstruction loss is the most natural
starting point for a single-run analysis.
[`get_best_run()`](https://gregorlueg.github.io/bixverse/reference/get_best_run.md)
returns it as an `NmfResult`, identical in shape to the output of
[`nmf_sc()`](https://gregorlueg.github.io/bixverse/reference/nmf_sc.md).

``` r

best_run <- get_best_run(t_cell_nmf_results_stabilised)
```

### Downstream analysis

What you do next depends entirely on what you are after. Some natural
directions:

- **Inspect the components**: rank genes by loading in `W` to read off
  what each component represents. Components that recover known
  programmes (cytotoxicity, exhaustion, ribosomal, IFN response, etc.)
  are a good sign; components that load on uninteresting genes you can
  ignore or merge.
- **Score cells**: `H` gives you per-cell activations for each
  programme. Project them onto your UMAP, correlate against existing
  annotations, or use them as features for downstream clustering.
- **Add activations to the obs table**:
  [`get_data()`](https://gregorlueg.github.io/bixverse/reference/get_data.md)
  on an `NmfResult` returns a data.table of cell activations that slots
  into the SingleCells obs table via
  [`add_sc_new_obs()`](https://gregorlueg.github.io/bixverse/reference/add_sc_new_obs.md),
  the same way module scores did earlier.
- **Consensus across runs**: with
  [`stabilised_nmf_sc()`](https://gregorlueg.github.io/bixverse/reference/stabilised_nmf_sc.md)
  you have `k * n_runs` W columns to play with. The classical move is to
  cluster these to identify components that recur across initialisations
  (the stable programmes) and components that appear only once (likely
  overfit to a particular seed).

We leave the choice of analysis here to you - the right downstream
depends on the biology and your questions to the data, not the API.

## Clean up

``` r

unlink(tempdir_pbmc, recursive = TRUE, force = TRUE)
```
