# Plot a consensus NMF k sweep

Stability and relative reconstruction error against `k`, the standard
cNMF diagnostic. Pick the `k` where stability is still high and the
error curve has not yet flattened out. Values of `k` whose consensus
step failed have no stability point.

## Usage

``` r
# S3 method for class 'NmfKSweepResult'
plot(x, ...)
```

## Arguments

- x:

  `NmfKSweepResult` object.

- ...:

  Additional params. Currently unused.

## Value

A `ggplot2` object with the two curves side by side.

## References

Kotliar et al., eLife, 2019

## Examples

``` r
# stability against reconstruction error across k
sc <- demo_single_cells()
res <- nmf_k_sweep_sc(
  sc,
  k_range = 2:4,
  n_runs = 3L,
  nmf_consensus_params = params_nmf_consensus(density_threshold = 2),
  .verbose = FALSE
)
plot(res)


unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
