# Wrapper function for consensus NMF parameters

Controls the clustering step of consensus NMF: the components of every
restart are pooled, outliers are dropped by local density, and the
survivors are k-means clustered into `k` groups whose median becomes the
consensus factor.

## Usage

``` r
params_nmf_consensus(
  consensus_target = c("h", "w"),
  n_neighbours = 0L,
  density_threshold = 0.5,
  kmeans_iters = 100L,
  kmeans_n_init = 3L
)
```

## Arguments

- consensus_target:

  String. One of `c("h", "w")`. `"h"` clusters the gene programmes (the
  spectra), which is what cNMF does and what you almost always want.
  `"w"` clusters in sample space instead, which on single cell data
  means cell space: it pools a dense `(k * n_runs) x n_cells` matrix and
  runs an exhaustive cosine search over it, so it gets expensive fast.

- n_neighbours:

  Integer. Neighbours used for the local density estimate. `0L` picks
  `ceiling(0.3 * n_runs)` for you.

- density_threshold:

  Numeric. Components whose mean cosine distance to their neighbours
  exceeds this are dropped as unstable. Cosine distance cannot exceed 2,
  so any value `>= 2` disables the filter entirely.

- kmeans_iters:

  Integer. Maximum k-means iterations.

- kmeans_n_init:

  Integer. Number of k-means restarts.

## Value

A list with the consensus NMF parameters.

## References

Kotliar et al., eLife, 2019
