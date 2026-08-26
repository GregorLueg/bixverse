# Sweep k for consensus NMF on a BulkCoExp

Runs the consensus clustering step across a range of `k` and reports
stability against reconstruction error, without keeping any factors.
Pick the `k` where stability is still high and the error curve has not
yet flattened out, then run
[`consensus_nmf_bulk()`](https://gregorlueg.github.io/bixverse/reference/consensus_nmf_bulk.md)
there.

## Usage

``` r
nmf_k_sweep_bulk(
  object,
  k_range,
  n_runs = 30L,
  preprocessing = c("none", "sd", "sqrt_sd"),
  nmf_hals_params = params_nmf_hals(),
  nmf_consensus_params = params_nmf_consensus(),
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

- k_range:

  Integer vector. The ranks to evaluate. Every entry must be at least 2.

- n_runs:

  Integer. Number of random restarts per `k`. Must be at least 2.

- preprocessing:

  String. One of `c("none", "sd", "sqrt_sd")`.

- nmf_hals_params:

  List. Output of
  [`params_nmf_hals()`](https://gregorlueg.github.io/bixverse/reference/params_nmf_hals.md).

- nmf_consensus_params:

  List. Output of
  [`params_nmf_consensus()`](https://gregorlueg.github.io/bixverse/reference/params_nmf_consensus.md).

- seed:

  Integer. Base random seed.

- .verbose:

  Boolean or integer `0L`/`1L`/`2L`. Controls verbosity.

## Value

An `NmfKSweepResult`, which is a data.table with one row per `k`.

## Details

This is a diagnostic, so it leaves the object alone and hands you the
result back directly.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on the returned
object gives you the two curves.

Cost is `length(k_range) * n_runs` full NMF fits, so keep both modest on
a first pass.

## References

Kotliar et al., eLife, 2019
