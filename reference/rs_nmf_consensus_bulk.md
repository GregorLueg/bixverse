# Run consensus NMF on a bulk expression matrix

**\[experimental\]** Runs `n_runs` HALS-NMF restarts, pools their
components, drops unstable ones by local density, k-means clusters the
survivors and refits the partner factor against the per-cluster median.
Expects samples x features.

## Usage

``` r
rs_nmf_consensus_bulk(
  x,
  k,
  preprocessing,
  nmf_hals_params,
  nmf_consensus_params,
  n_runs,
  seed,
  verbose
)
```

## Arguments

- x:

  Numerical matrix. Rows = samples, columns = features.

- k:

  Integer. Number of latent factors. Must be at least 2.

- preprocessing:

  String. One of `c("none", "sd", "sqrt_sd")`.

- nmf_hals_params:

  Named list. See
  [`params_nmf_hals()`](https://gregorlueg.github.io/bixverse/reference/params_nmf_hals.md).
  The `nmf_init` field is ignored, restarts always use random
  initialisation.

- nmf_consensus_params:

  Named list. See
  [`params_nmf_consensus()`](https://gregorlueg.github.io/bixverse/reference/params_nmf_consensus.md).

- n_runs:

  Integer. Number of restarts. Must be at least 2.

- seed:

  Integer. Base random seed. Restart `i` uses `seed + i`.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following items

- w - The `W` matrix of shape `n_samples x k`.

- h - The `H` matrix of shape `k x n_features`.

- rel_error - Reconstruction error relative to the squared Frobenius
  norm of the input. Not comparable with the absolute `final_loss` the
  single-run version returns.

- rel_run_errors - The same, per restart.

- labels - Integer vector of length `k * n_runs`. Cluster each pooled
  component landed in, `NA` if it was dropped.

- local_density - Mean cosine distance to the nearest neighbours per
  pooled component. Plot this to pick a `density_threshold`.

- kept - 1-indexed positions of the surviving pooled components.

- silhouette - Silhouette per survivor, aligned with `kept`.

- stability - Mean silhouette over the survivors.

- cluster_sizes - Number of survivors per cluster.

- n_dropped - Number of pooled components removed.

- n_empty_clusters - Number of clusters left with no members.

## References

Kotliar et al., eLife, 2019
