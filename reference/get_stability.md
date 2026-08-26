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
