# Constructor for stabilised (multi-run) NMF results

Stores the column-bound W matrices across runs, per-run H matrices,
per-run losses and convergence flags, plus the index of the best run.

## Usage

``` r
new_stabilised_nmf_result(
  nmf_res,
  gene_ids,
  cell_ids,
  cell_indices,
  source_class,
  params
)
```

## Arguments

- nmf_res:

  Named list. Output of
  [`rs_nmf_multi_sc()`](https://gregorlueg.github.io/bixverse/reference/rs_nmf_multi_sc.md)
  or
  [`rs_nmf_multi_mc()`](https://gregorlueg.github.io/bixverse/reference/rs_nmf_multi_mc.md).
  Must contain `w_all`, `h_per_run`, `losses`, `converged`, `best_idx`.

- gene_ids:

  Character vector. Gene identifiers for the rows of W.

- cell_ids:

  Character vector. Cell identifiers for the columns of H.

- cell_indices:

  Integer vector. 0-indexed cell positions used in the run.

- source_class:

  String. One of `c("SingleCells", "MetaCells")`.

- params:

  List. The full set of parameters used for the run.

## Value

An object of class `StabilisedNmfResult`.
