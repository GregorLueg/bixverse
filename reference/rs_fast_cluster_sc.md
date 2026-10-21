# Runs fast clustering on the data

**\[experimental\]** Coarsens the cells with k-means, builds a graph
over the centroids (a kNN graph with one of the symmetrisations, or the
cell-mediated graph) and runs Louvain or Leiden on it. The community of
each centroid is propagated back to its cells. With `n_seeds > 1` every
resolution is run once per seed, the labels come from the seed with the
lowest mean conductance and stability metrics across the seeds are
returned.

## Usage

``` r
rs_fast_cluster_sc(embd, resolutions, fc_params, n_seeds, seed, verbose)
```

## Arguments

- embd:

  Numeric matrix. The embedding, cells x features.

- resolutions:

  Numeric vector. The resolutions to iterate through.

- fc_params:

  Named list. The fast clustering parameters, see
  [`params_sc_fast_cluster()`](https://gregorlueg.github.io/bixverse/reference/params_sc_fast_cluster.md).
  Also carries the k-means keys (`km_type`, `kmeans_iters`,
  `batch_size`, `drift_threshold`, `lr_alpha`).

- n_seeds:

  Integer. Number of seeds per resolution. `1L` runs a single pass
  without stability metrics.

- seed:

  Integer. For reproducibility.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following elements:

- membership - List with one integer membership vector per resolution
  (0-indexed).

- k_means_cluster - Integer vector with the k-means centroid per cell
  (0-indexed).

- centroids - Numerical matrix with the k-means centroids.

- stats - `NULL` for `n_seeds = 1L`, otherwise a list with `mean_ari`,
  `median_ari`, `mean_conductance`, `median_conductance` and
  `mean_n_comms`, one value per resolution.
