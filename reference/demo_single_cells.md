# Ready-made `SingleCells` object for examples and tests

Wires
[`generate_single_cell_test_data()`](https://gregorlueg.github.io/bixverse/reference/generate_single_cell_test_data.md)
into a `SingleCells` object on disk in one call, so examples do not have
to repeat the whole ingestion dance. The default is deliberately tiny
(500 cells x 50 genes) and the quality thresholds are loose enough that
every cell survives. This is synthetic data for demonstration and
testing, not something to analyse.

## Usage

``` r
demo_single_cells(
  dir = tempfile("bixverse_demo"),
  prepped = TRUE,
  syn_data_params = params_sc_synthetic_data(n_cells = 500L, n_genes = 50L),
  hvg_no = 30L,
  no_pcs = 10L,
  k = 15L,
  seed = 42L,
  .verbose = FALSE
)
```

## Arguments

- dir:

  String. Directory to hold the object. Created if it does not exist.
  Defaults to a fresh path under the session
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html). Remove it with
  `unlink(dir, recursive = TRUE)` when you are done.

- prepped:

  Boolean. Run the standard HVG -\> PCA -\> kNN chain before returning?
  Defaults to `TRUE`.

- syn_data_params:

  List. Parameters for the synthetic data, see
  [`params_sc_synthetic_data()`](https://gregorlueg.github.io/bixverse/reference/params_sc_synthetic_data.md).

- hvg_no:

  Integer. Number of highly variable genes, `prepped` only.

- no_pcs:

  Integer. Number of principal components, `prepped` only.

- k:

  Integer. Number of nearest neighbours, `prepped` only.

- seed:

  Integer. Seed for the data generation.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

A `SingleCells` object backed by `dir`.

## Examples

``` r
# a prepped object, ready for clustering
sc <- demo_single_cells()
sc <- find_clusters_sc(sc, res = 1.0)
sc
#> Single cell experiment (Single Cells).
#>   No cells (original): 500
#>    To keep n: 500
#>   No genes: 50
#>   HVG calculated: TRUE
#>   PCA calculated: TRUE
#>   Other embeddings: none
#>   KNN generated: TRUE
#>   SNN generated: TRUE
#>   MAGIC imputed: none
#>   Residual model: none
#>   Stale artefacts: none

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
