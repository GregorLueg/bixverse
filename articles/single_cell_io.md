# Single-cell data I/O in bixverse

## Intro

`bixverse` keeps single-cell counts on disk – in compact Rust binaries
for the count matrix and a DuckDB for the `obs` / `var` metadata – so
that millions of cells can be analysed on a laptop. Ingestion happens
once, up front, through format-specific readers that stream the data
into that on-disk layout and apply the basic cell/gene QC and
normalisations along the way.

This vignette walks through every reader on a small synthetic data set,
so it runs end-to-end without any downloads. It covers:

- in-memory R objects (`load_r_data`)
- 10x CellRanger MTX, single and multi-sample (`load_mtx`,
  `load_multi_mtx`)
- AnnData `h5ad`, single and multi-sample (`load_h5ad`,
  `load_multi_h5ad`)
- 10x CellRanger HDF5, single and multi-sample (`load_tenx_h5`,
  `load_multi_tenx_h5`)
- checkpointing to and from disk (`save_sc_exp_to_disk`,
  `load_existing`)

If you have not seen the `SingleCells` class before, read the [design
choices](https://gregorlueg.github.io/bixverse/articles/design_single_cell.html)
vignette first.

``` r

library(bixverse)
library(data.table)
#> 
#> Attaching package: 'data.table'
#> The following object is masked from 'package:base':
#> 
#>     %notin%
library(Matrix)
library(magrittr)
```

## A small synthetic data set

We generate a synthetic data set and treat its two halves as two
“samples”, so that the multi-sample readers can be demonstrated too. The
counts are a cell-by-gene `dgRMatrix`; `obs` and `var` are
`data.table`s.

``` r

sim <- generate_single_cell_test_data()

counts <- sim$counts # cells x genes (dgRMatrix)
obs <- sim$obs
var <- sim$var
dim(counts)
#> [1] 1000  100

# two "samples" for the multi-sample examples
idx_a <- 1:500
idx_b <- 501:1000
counts_a <- as(counts[idx_a, ], "RsparseMatrix")
counts_b <- as(counts[idx_b, ], "RsparseMatrix")

# scratch directories (one on-disk store per loaded object)
base <- file.path(tempdir(), "bixverse_io")
stores <- c(
  "mtx_a",
  "mtx_b",
  "mtx_swapped",
  "store_r",
  "store_mtx",
  "store_mtx_swapped",
  "store_multi_mtx",
  "store_h5",
  "store_multi_h5",
  "store_10x",
  "store_multi_10x"
)
for (s in stores) {
  dir.create(file.path(base, s), recursive = TRUE, showWarnings = FALSE)
}

# light QC thresholds, suitable for this small toy data set
qc <- params_sc_min_quality(
  min_unique_genes = 5L,
  min_lib_size = 10L,
  min_cells = 3L
)
```

## In-memory R objects

The simplest entry point: hand the reader a cell-by-gene `dgRMatrix`
plus the `obs` and `var` tables directly.

``` r

sc_r <- SingleCells(dir_data = file.path(base, "store_r"))
sc_r <- load_r_data(
  object = sc_r,
  counts = counts,
  obs = obs,
  var = var,
  sc_qc_param = qc,
  .verbose = FALSE
)
sc_r
#> Single cell experiment (Single Cells).
#>   No cells (original): 1000
#>    To keep n: 1000
#>   No genes: 100
#>   HVG calculated: FALSE
#>   PCA calculated: FALSE
#>   Other embeddings: none
#>   KNN generated: FALSE
#>   SNN generated: FALSE
```

## 10x CellRanger MTX

### Single sample

We write one sample out in CellRanger flat-file form and point an
explicit
[`params_sc_mtx_io()`](https://gregorlueg.github.io/bixverse/reference/params_sc_mtx_io.md)
at the trio (for real CellRanger output, the convenience helper
[`get_cell_ranger_params()`](https://gregorlueg.github.io/bixverse/reference/get_cell_ranger_params.md)
auto-locates `matrix.mtx` / `barcodes.tsv` / `genes.tsv`).

Two
[`params_sc_mtx_io()`](https://gregorlueg.github.io/bixverse/reference/params_sc_mtx_io.md)
arguments handle the common deviations from the strict CellRanger
convention:

- `cells_as_rows` – whether the .mtx has cells on the rows. The
  CellRanger default is `FALSE` (genes are rows). Set to `TRUE` for
  externally produced matrices that are transposed.
- `has_hdr` – whether the barcodes/features flat files carry a header
  row. CellRanger itself writes headerless TSVs (`has_hdr = FALSE`), but
  CSV/TSV exports from downstream pipelines often carry one. The example
  below uses `has_hdr = TRUE` because
  [`write_cellranger_output()`](https://gregorlueg.github.io/bixverse/reference/write_cellranger_output.md)
  emits headers.

``` r

dir_a <- file.path(base, "mtx_a")
# real CellRanger MTX is genes x cells, matching cells_as_rows = FALSE below
write_cellranger_output(
  dir_a,
  counts_a,
  obs[idx_a],
  var,
  rows = "genes",
  .verbose = FALSE
)

mtx_params <- params_sc_mtx_io(
  path_mtx = file.path(dir_a, "mat.mtx"),
  path_obs = file.path(dir_a, "barcodes.csv"),
  path_var = file.path(dir_a, "features.csv"),
  cells_as_rows = FALSE,
  has_hdr = TRUE
)

sc_mtx <- SingleCells(dir_data = file.path(base, "store_mtx"))
sc_mtx <- load_mtx(
  object = sc_mtx,
  sc_mtx_io_param = mtx_params,
  sc_qc_param = qc,
  mtx_streaming = FALSE,
  .verbose = FALSE
)
sc_mtx
#> Single cell experiment (Single Cells).
#>   No cells (original): 500
#>    To keep n: 500
#>   No genes: 100
#>   HVG calculated: FALSE
#>   PCA calculated: FALSE
#>   Other embeddings: none
#>   KNN generated: FALSE
#>   SNN generated: FALSE
```

### Multiple samples

[`prescan_mtx_dirs()`](https://gregorlueg.github.io/bixverse/reference/prescan_mtx_dirs.md)
takes one directory per sample, builds the shared gene universe across
them, and
[`load_multi_mtx()`](https://gregorlueg.github.io/bixverse/reference/load_multi_mtx.md)
then stacks all cells into a single object with a global gene QC. The
originating sample is recorded in the `exp_id` `obs` column.

``` r

dir_b <- file.path(base, "mtx_b")
write_cellranger_output(
  dir_b,
  counts_b,
  obs[idx_b],
  var,
  rows = "genes",
  .verbose = FALSE
)

ps_mtx <- prescan_mtx_dirs(
  dirs = c(dir_a, dir_b),
  exp_ids = c("sample_a", "sample_b"),
  has_hdr = TRUE
)

sc_multi_mtx <- SingleCells(dir_data = file.path(base, "store_multi_mtx"))
sc_multi_mtx <- load_multi_mtx(
  object = sc_multi_mtx,
  prescan_result = ps_mtx,
  sc_qc_param = qc,
  .verbose = FALSE
)
sc_multi_mtx
#> Single cell experiment (Single Cells).
#>   No cells (original): 1000
#>    To keep n: 1000
#>   No genes: 100
#>   HVG calculated: FALSE
#>   PCA calculated: FALSE
#>   Other embeddings: none
#>   KNN generated: FALSE
#>   SNN generated: FALSE

# cells per sample
get_sc_obs(sc_multi_mtx)[, .N, by = exp_id]
#>      exp_id     N
#>      <char> <int>
#> 1: sample_a   500
#> 2: sample_b   500
```

### File layout variations

To make the `cells_as_rows` and `has_hdr` knobs concrete, here is the
same sample written with cells on the rows of the .mtx (i.e. transposed
relative to CellRanger). The only thing that changes at load time is
`cells_as_rows = TRUE`:

``` r

dir_swapped <- file.path(base, "mtx_swapped")
write_cellranger_output(
  dir_swapped,
  counts_a,
  obs[idx_a],
  var,
  rows = "cells",
  .verbose = FALSE
)

mtx_params_swapped <- params_sc_mtx_io(
  path_mtx = file.path(dir_swapped, "mat.mtx"),
  path_obs = file.path(dir_swapped, "barcodes.csv"),
  path_var = file.path(dir_swapped, "features.csv"),
  cells_as_rows = TRUE,
  has_hdr = TRUE
)

sc_mtx_swapped <- SingleCells(dir_data = file.path(base, "store_mtx_swapped"))
sc_mtx_swapped <- load_mtx(
  object = sc_mtx_swapped,
  sc_mtx_io_param = mtx_params_swapped,
  sc_qc_param = qc,
  mtx_streaming = FALSE,
  .verbose = FALSE
)
sc_mtx_swapped
#> Single cell experiment (Single Cells).
#>   No cells (original): 500
#>    To keep n: 500
#>   No genes: 100
#>   HVG calculated: FALSE
#>   PCA calculated: FALSE
#>   Other embeddings: none
#>   KNN generated: FALSE
#>   SNN generated: FALSE
```

The resulting on-disk object is identical to the standard-layout
single-sample load above; only the path through the reader differs.

For headerless barcodes/features files (the strict CellRanger
convention), set `has_hdr = FALSE`. The same flag applies uniformly
across all input directories in
[`prescan_mtx_dirs()`](https://gregorlueg.github.io/bixverse/reference/prescan_mtx_dirs.md)
for the multi-sample case.

## AnnData h5ad

### Single

[`read_h5ad_metadata()`](https://gregorlueg.github.io/bixverse/reference/read_h5ad_metadata.md)
lets you inspect the `obs` / `var` and dimensions without loading the
counts;
[`load_h5ad()`](https://gregorlueg.github.io/bixverse/reference/load_h5ad.md)
ingests the matrix.

``` r

h5_a <- file.path(base, "sample_a.h5ad")
write_h5ad_sc(h5_a, counts_a, obs[idx_a], var, .verbose = FALSE)

meta <- read_h5ad_metadata(h5_a)
meta$dims
#> obs var 
#> 500 100

sc_h5 <- SingleCells(dir_data = file.path(base, "store_h5"))
sc_h5 <- load_h5ad(
  object = sc_h5,
  h5_path = h5_a,
  sc_qc_param = qc,
  .verbose = FALSE
)
sc_h5
#> Single cell experiment (Single Cells).
#>   No cells (original): 500
#>    To keep n: 500
#>   No genes: 100
#>   HVG calculated: FALSE
#>   PCA calculated: FALSE
#>   Other embeddings: none
#>   KNN generated: FALSE
#>   SNN generated: FALSE
```

### Multiple

[`prescan_h5ad_files()`](https://gregorlueg.github.io/bixverse/reference/prescan_h5ad_files.md)
accepts a named vector of paths (the names become the sample
identifiers);
[`load_multi_h5ad()`](https://gregorlueg.github.io/bixverse/reference/load_multi_h5ad.md)
mirrors the multi-MTX workflow.

``` r

h5_b <- file.path(base, "sample_b.h5ad")
write_h5ad_sc(h5_b, counts_b, obs[idx_b], var, .verbose = FALSE)

ps_h5 <- prescan_h5ad_files(c(sample_a = h5_a, sample_b = h5_b))

sc_multi_h5 <- SingleCells(dir_data = file.path(base, "store_multi_h5"))
sc_multi_h5 <- load_multi_h5ad(
  object = sc_multi_h5,
  prescan_result = ps_h5,
  sc_qc_param = qc,
  .verbose = FALSE
)
sc_multi_h5
#> Single cell experiment (Single Cells).
#>   No cells (original): 1000
#>    To keep n: 1000
#>   No genes: 100
#>   HVG calculated: FALSE
#>   PCA calculated: FALSE
#>   Other embeddings: none
#>   KNN generated: FALSE
#>   SNN generated: FALSE
```

## 10x CellRanger HDF5

### Single sample

[`write_tenx_h5_sc()`](https://gregorlueg.github.io/bixverse/reference/write_tenx_h5_sc.md)
expects the counts plus a barcode vector and a `features` table (with at
least `id` and `name`).
[`load_tenx_h5()`](https://gregorlueg.github.io/bixverse/reference/load_tenx_h5.md)
reads the gene-expression modality (V3 files may also carry
e.g. Antibody Capture, which is filtered out via `feature_type`).

``` r

features <- data.table(id = var$gene_id, name = var$ensembl_id)
h5_tenx <- file.path(base, "sample_a_10x.h5")
write_tenx_h5_sc(
  f_path = h5_tenx,
  counts = counts_a,
  barcodes = obs[idx_a]$cell_id,
  features = features,
  version = "v3"
)

sc_10x <- SingleCells(dir_data = file.path(base, "store_10x"))
sc_10x <- load_tenx_h5(
  object = sc_10x,
  h5_path = h5_tenx,
  sc_qc_param = qc,
  feature_type = "Gene Expression",
  .verbose = FALSE
)
sc_10x
#> Single cell experiment (Single Cells).
#>   No cells (original): 500
#>    To keep n: 500
#>   No genes: 100
#>   HVG calculated: FALSE
#>   PCA calculated: FALSE
#>   Other embeddings: none
#>   KNN generated: FALSE
#>   SNN generated: FALSE
```

### Multiple samples

[`prescan_tenx_h5_files()`](https://gregorlueg.github.io/bixverse/reference/prescan_tenx_h5_files.md)
takes a named vector of `.h5` paths (the names become the sample
identifiers), inspects each file for its CellRanger version (v2 vs v3),
filters to the requested `feature_type`, and builds the global gene
space.
[`load_multi_tenx_h5()`](https://gregorlueg.github.io/bixverse/reference/load_multi_tenx_h5.md)
then mirrors the multi-MTX / multi-h5ad workflow.

`gene_universe` controls how the feature space is reconciled across
files:

- `"intersection"` keeps genes present in **every** input.
- `"union"` keeps genes present in **any** input; cells from a file that
  lacks a given gene get zero counts for it.

Mixed v2 / v3 inputs in the same call are handled transparently, and
non-gene modalities (e.g. Antibody Capture) are silently dropped when
`feature_type = "Gene Expression"`.

``` r

features_b <- data.table(
  id = var$gene_id,
  name = var$ensembl_id,
  feature_type = "Gene Expression"
)
h5_tenx_b <- file.path(base, "sample_b_10x.h5")
write_tenx_h5_sc(
  f_path = h5_tenx_b,
  counts = counts_b,
  barcodes = obs[idx_b]$cell_id,
  features = features_b,
  version = "v3"
)

ps_10x <- prescan_tenx_h5_files(
  h5_paths = c(sample_a = h5_tenx, sample_b = h5_tenx_b),
  feature_type = "Gene Expression",
  gene_universe = "intersection",
  .verbose = FALSE
)

sc_multi_10x <- SingleCells(dir_data = file.path(base, "store_multi_10x"))
sc_multi_10x <- load_multi_tenx_h5(
  object = sc_multi_10x,
  prescan_result = ps_10x,
  sc_qc_param = qc,
  .verbose = FALSE
)
sc_multi_10x
#> Single cell experiment (Single Cells).
#>   No cells (original): 1000
#>    To keep n: 1000
#>   No genes: 100
#>   HVG calculated: FALSE
#>   PCA calculated: FALSE
#>   Other embeddings: none
#>   KNN generated: FALSE
#>   SNN generated: FALSE

# cells per sample
get_sc_obs(sc_multi_10x)[, .N, by = exp_id]
#>      exp_id     N
#>      <char> <int>
#> 1: sample_a   500
#> 2: sample_b   500
```

For multi-modal `.h5` files that carry both Gene Expression and Antibody
Capture features, load the gene side as above and read the ADT side
separately with
[`read_tenx_h5_adt()`](https://gregorlueg.github.io/bixverse/reference/read_tenx_h5_adt.md),
then attach it to a `SingleCellsMultiModal` object via
[`add_adt_counts_sc()`](https://gregorlueg.github.io/bixverse/reference/add_adt_counts_sc.md).

## Checkpointing to disk

Once an object is loaded, the in-memory mappings can be persisted so a
fresh session can be restored without re-ingesting the counts. The
on-disk binary and DuckDB already live in `dir_data`;
[`save_sc_exp_to_disk()`](https://gregorlueg.github.io/bixverse/reference/save_sc_exp_to_disk.md)
adds the lightweight mapping cache.

``` r

save_sc_exp_to_disk(sc_h5, type = "rds")

sc_restored <- SingleCells(dir_data = file.path(base, "store_h5"))
sc_restored <- load_existing(sc_restored)
#> Found stored data from save_sc_exp_to_disk(). Loading that one into the object.
sc_restored
#> Single cell experiment (Single Cells).
#>   No cells (original): 500
#>    To keep n: 500
#>   No genes: 100
#>   HVG calculated: FALSE
#>   PCA calculated: FALSE
#>   Other embeddings: none
#>   KNN generated: FALSE
#>   SNN generated: FALSE
```

## Summary

| Format | Single sample | Multiple samples |
|----|----|----|
| In-memory R | [`load_r_data()`](https://gregorlueg.github.io/bixverse/reference/load_r_data.md) | – |
| 10x MTX | [`load_mtx()`](https://gregorlueg.github.io/bixverse/reference/load_mtx.md) | [`prescan_mtx_dirs()`](https://gregorlueg.github.io/bixverse/reference/prescan_mtx_dirs.md) + [`load_multi_mtx()`](https://gregorlueg.github.io/bixverse/reference/load_multi_mtx.md) |
| AnnData h5ad | [`load_h5ad()`](https://gregorlueg.github.io/bixverse/reference/load_h5ad.md) | [`prescan_h5ad_files()`](https://gregorlueg.github.io/bixverse/reference/prescan_h5ad_files.md) + [`load_multi_h5ad()`](https://gregorlueg.github.io/bixverse/reference/load_multi_h5ad.md) |
| 10x HDF5 | [`load_tenx_h5()`](https://gregorlueg.github.io/bixverse/reference/load_tenx_h5.md) | [`prescan_tenx_h5_files()`](https://gregorlueg.github.io/bixverse/reference/prescan_tenx_h5_files.md) + [`load_multi_tenx_h5()`](https://gregorlueg.github.io/bixverse/reference/load_multi_tenx_h5.md) |
| Seurat | [`load_seurat()`](https://gregorlueg.github.io/bixverse/reference/load_seurat.md) | – |

All readers share the same `sc_qc_param` QC interface and write into the
same on-disk `SingleCells` layout, so everything downstream
(normalisation, HVGs, PCA, neighbours, clustering, …) is identical
regardless of how the data came in.
