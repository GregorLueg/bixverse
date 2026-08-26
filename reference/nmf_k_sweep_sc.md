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
