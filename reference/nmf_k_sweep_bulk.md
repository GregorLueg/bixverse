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

## Examples

``` r
# small sweep, keep both the grid and the restarts modest
syn <- generate_gene_module_data(n_samples = 24L, n_genes = 60L)
# NMF needs a non-negative matrix
mat <- syn$data - min(syn$data)
obj <- BulkCoExp(mat, syn$meta_data)
obj <- preprocess_bulk_coexp(
  obj, hvg = NULL, scaling = FALSE, .verbose = FALSE
)
sweep_res <- nmf_k_sweep_bulk(
  obj, k_range = 3:5, n_runs = 5L, .verbose = FALSE
)
sweep_res
#> NmfKSweepResult (consensus NMF k sweep)
#>   Source class:     BulkCoExp
#>   k range:          3 to 5
#>   No runs per k:    5
#>   Most stable k:    3 (stability = 1)
#> 
#>        k stability  best_error median_error consensus_failed n_dropped
#>    <int>     <num>       <num>        <num>           <lgcl>     <int>
#> 1:     3 0.9999921 0.110235887  0.110236112            FALSE         0
#> 2:     4 0.9998356 0.003672419  0.003675283            FALSE         0
#> 3:     5 0.9234360 0.003480078  0.003504362            FALSE         0
#>    n_empty_clusters n_converged
#>               <int>       <int>
#> 1:                0           5
#> 2:                0           5
#> 3:                0           5
```
