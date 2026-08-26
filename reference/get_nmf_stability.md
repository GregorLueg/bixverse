# Get the multi-run NMF diagnostics

Getter function for the diagnostics of a multi-run NMF fit. For a
[`stabilised_nmf_bulk()`](https://gregorlueg.github.io/bixverse/reference/stabilised_nmf_bulk.md)
fit that is the per-run losses, convergence flags and the best-run
index. For a
[`consensus_nmf_bulk()`](https://gregorlueg.github.io/bixverse/reference/consensus_nmf_bulk.md)
fit it is the cluster stability, the relative errors and the
per-component clustering table. Returns `NULL` for a single-run fit.

## Usage

``` r
get_nmf_stability(object)
```

## Arguments

- object:

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

## Value

For a stabilised fit, a list with `losses`, `converged` and `best_idx`.
For a consensus fit, a list with `stability`, `rel_error`,
`rel_run_errors`, `clusters`, `cluster_sizes`, `n_dropped` and
`n_empty_clusters`. `NULL` if neither was run.
