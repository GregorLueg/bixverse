# Sweep k for consensus NMF on single cell or meta cell data

Runs the consensus clustering step across a range of `k` and reports
stability against reconstruction error, without keeping any factors.
Pick the `k` where stability is still high and the error curve has not
yet flattened out, then run
[`consensus_nmf_sc()`](https://gregorlueg.github.io/bixverse/reference/consensus_nmf_sc.md)
there.

## Usage

``` r
nmf_k_sweep_sc(
  object,
  k_range,
  cell_ids = NULL,
  gene_ids = NULL,
  preprocessing = "none",
  use_second_layer = TRUE,
  nmf_hals_params = params_nmf_hals(),
  nmf_consensus_params = params_nmf_consensus(),
  n_runs = 30L,
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` or `MetaCells` class.

- k_range:

  Integer vector. The ranks to evaluate. Every entry must be at least 2.

- cell_ids:

  Optional character. Cell ids (or meta cell ids) to restrict the NMF
  to. If `NULL`, uses
  [`get_cells_to_keep()`](https://gregorlueg.github.io/bixverse/reference/get_cells_to_keep.md)
  for `SingleCells` and all meta cells for `MetaCells`.

- gene_ids:

  Optional character. Gene ids to restrict the NMF to. If `NULL`, uses
  [`get_hvg()`](https://gregorlueg.github.io/bixverse/reference/get_hvg.md)
  on the object.

- preprocessing:

  String. One of `c("none", "sd", "sqrt_sd")`.

- use_second_layer:

  Boolean. If `TRUE`, runs NMF on the normalised counts (recommended);
  if `FALSE`, on the raw counts.

- nmf_hals_params:

  List, see
  [`params_nmf_hals()`](https://gregorlueg.github.io/bixverse/reference/params_nmf_hals.md).
  The `nmf_init` field is ignored, restarts always use random
  initialisation.

- nmf_consensus_params:

  List, see
  [`params_nmf_consensus()`](https://gregorlueg.github.io/bixverse/reference/params_nmf_consensus.md).

- n_runs:

  Integer. Number of random restarts. Must be at least 2.

- seed:

  Integer. Base random seed. Restart `i` uses `seed + i`, and the
  k-means step is seeded from it too.

- .verbose:

  Boolean or integer. Verbosity.

## Value

An `NmfKSweepResult`, which is a data.table with one row per `k`.

## Details

This is a diagnostic, so it leaves the object alone and hands the result
back directly. [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
on it gives you the two curves.

Cost is `length(k_range) * n_runs` full NMF fits. On the `SingleCells`
path the counts are read from disk once and reused across every `k`, but
the fits themselves are not free, so keep both modest on a first pass.

## References

Kotliar et al., eLife, 2019

## Examples

``` r
# stability against reconstruction error across three ranks
sc <- demo_single_cells()
res <- nmf_k_sweep_sc(
  sc,
  k_range = 2:4,
  n_runs = 5L,
  nmf_consensus_params = params_nmf_consensus(density_threshold = 2),
  .verbose = FALSE
)
res
#> NmfKSweepResult (consensus NMF k sweep)
#>   Source class:     SingleCells
#>   k range:          2 to 4
#>   No runs per k:    5
#>   Most stable k:    3 (stability = 0.9944)
#> 
#>        k stability best_error median_error consensus_failed n_dropped
#>    <int>     <num>      <num>        <num>           <lgcl>     <int>
#> 1:     2 0.9943518  0.2972628    0.2973181            FALSE         0
#> 2:     3 0.9944324  0.2518355    0.2518508            FALSE         0
#> 3:     4 0.9003837  0.2320080    0.2320648            FALSE         0
#>    n_empty_clusters n_converged
#>               <int>       <int>
#> 1:                0           5
#> 2:                0           5
#> 3:                0           5

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
