# Get the consensus NMF stability diagnostics

Returns the mean silhouette of the consensus clusters, the relative
reconstruction errors and the per-component clustering table. A
stability near 1 means every restart found the same programmes; a low
one means the factorisation is not reproducible at this `k`.

## Usage

``` r
get_stability(x)

# S3 method for class 'ConsensusNmfResult'
get_stability(x)
```

## Arguments

- x:

  An object holding consensus NMF results.

## Value

A list with `stability`, `rel_error`, `rel_run_errors`, `clusters`,
`cluster_sizes`, `n_dropped` and `n_empty_clusters`.

## Examples

``` r
# how much the restarts agreed on the consensus factors
sc <- demo_single_cells()
res <- consensus_nmf_sc(
  sc,
  k = 5L,
  n_runs = 5L,
  nmf_consensus_params = params_nmf_consensus(density_threshold = 2),
  .verbose = FALSE
)
get_stability(res)$stability
#> [1] 0.9446008

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
