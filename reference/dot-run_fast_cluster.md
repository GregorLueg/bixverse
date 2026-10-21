# Run the fast clustering on an embedding and wrap the results

Shared by the `SingleCells` and `SingleCellsSubset` methods of
[`fast_cluster_sc()`](https://gregorlueg.github.io/bixverse/reference/fast_cluster_sc.md).
Resolves the `NULL` defaults of
[`params_sc_fast_cluster()`](https://gregorlueg.github.io/bixverse/reference/params_sc_fast_cluster.md)
before handing over to Rust.

## Usage

``` r
.run_fast_cluster(
  embd,
  cells_to_use,
  resolutions,
  fc_params,
  n_seeds,
  return_kmeans,
  seed,
  .verbose
)
```

## Arguments

- embd:

  Numeric matrix. Cells x features.

- cells_to_use:

  Integer vector. 0-indexed cell indices of the rows of `embd`.

- resolutions:

  Numeric vector. Modularity resolutions.

- fc_params:

  List. Output of
  [`params_sc_fast_cluster()`](https://gregorlueg.github.io/bixverse/reference/params_sc_fast_cluster.md).

- n_seeds:

  Integer. Number of seeds per resolution. `1L` runs a single pass
  without stability statistics.

- return_kmeans:

  Boolean. Return k-means assignments and centroids.

- seed:

  Integer. Reproducibility.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

`SingleCellFastClusters` S3 object, see
[`fast_cluster_sc()`](https://gregorlueg.github.io/bixverse/reference/fast_cluster_sc.md).
