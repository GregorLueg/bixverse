# Build a Symphony reference from a SingleCells object

Runs sparse PCA over the provided HVGs, runs Harmony for batch
correction, and compresses the result into the cached terms used at
query time. The Harmony version is auto-detected from the class of
`harmony_params`. Optionally snapshots one or more obs columns (e.g.
cell type annotations) into the reference's `labels` slot for downstream
label transfer.

## Usage

``` r
build_symphony_ref(
  object,
  batch_column,
  additional_batch_columns = NULL,
  hvg,
  harmony_params = params_sc_harmony(),
  pca_params = params_sc_pca(),
  no_pcs = 30L,
  slim = FALSE,
  label_columns = NULL,
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` (the reference).

- batch_column:

  String. Primary batch column in the obs table.

- additional_batch_columns:

  Optional character vector.

- hvg:

  Integer vector. R-style 1-based HVG indices into the gene universe of
  `object`. Must be provided explicitly.

- harmony_params:

  List. Output of
  [`params_sc_harmony()`](https://gregorlueg.github.io/bixverse/reference/params_sc_harmony.md)
  or
  [`params_sc_harmony_v2()`](https://gregorlueg.github.io/bixverse/reference/params_sc_harmony_v2.md).

- pca_params:

  List. Output of
  [`params_sc_pca()`](https://gregorlueg.github.io/bixverse/reference/params_sc_pca.md).

- no_pcs:

  Integer. Number of principal components.

- slim:

  Boolean. If `TRUE`, drops `z_orig` and `r` from the returned
  reference. `z_corr` is always kept (needed for kNN label transfer).

- label_columns:

  Optional character vector of obs column names to snapshot into the
  reference. Read in `cells_to_keep` order to align with `z_corr`.
  Required for
  [`transfer_labels_symphony()`](https://gregorlueg.github.io/bixverse/reference/transfer_labels_symphony.md)
  unless populated later via
  [`add_symphony_labels()`](https://gregorlueg.github.io/bixverse/reference/add_symphony_labels.md).

- seed:

  Integer.

- .verbose:

  Boolean or integer.

## Value

A
[SymphonyReference](https://gregorlueg.github.io/bixverse/reference/SymphonyReference.md)
object.

## Examples

``` r
# two batch reference, cell type labels snapshotted along the way
ref <- demo_single_cells(
  syn_data_params = params_sc_synthetic_data(
    n_cells = 500L,
    n_genes = 50L,
    n_batches = 2L,
    batch_effect_strength = "medium"
  )
)
symphony_ref <- build_symphony_ref(
  ref,
  batch_column = "batch_index",
  hvg = get_hvg(ref) + 1L,
  harmony_params = params_sc_harmony(k = 10L),
  no_pcs = 10L,
  label_columns = "cell_grp",
  .verbose = FALSE
)
symphony_ref
#> Symphony reference
#>   Harmony backend: v1
#>   No HVGs: 30
#>   No PCs: 10
#>   No clusters: 10
#>   Batch variables: batch_index
#>   Slim: FALSE
#>   Labels: cell_grp

unlink(ref@dir_data, recursive = TRUE, force = TRUE)
```
