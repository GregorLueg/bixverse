# Constructor for NMF k sweep results

Wraps the diagnostics of a consensus NMF sweep over several `k` in a
data.table, one row per `k`. Subclassing the data.table rather than
boxing it means `[`, [`order()`](https://rdrr.io/r/base/order.html) and
friends keep working.

## Usage

``` r
new_nmf_k_sweep_result(sweep_res, source_class, params)
```

## Arguments

- sweep_res:

  Named list. Output of
  [`rs_nmf_k_sweep_bulk()`](https://gregorlueg.github.io/bixverse/reference/rs_nmf_k_sweep_bulk.md),
  [`rs_nmf_k_sweep_sc()`](https://gregorlueg.github.io/bixverse/reference/rs_nmf_k_sweep_sc.md)
  or
  [`rs_nmf_k_sweep_mc()`](https://gregorlueg.github.io/bixverse/reference/rs_nmf_k_sweep_mc.md).
  Must contain `k`, `stability`, `best_error`, `median_error`,
  `consensus_failed`, `n_dropped`, `n_empty_clusters` and `n_converged`.

- source_class:

  String. One of `c("BulkCoExp", "SingleCells", "MetaCells")`.

- params:

  List. The full set of parameters used for the sweep.

## Value

An object of class `NmfKSweepResult`, which is also a data.table.
`stability` is `NA` for any `k` where the consensus step failed.

## References

Kotliar et al., eLife, 2019
