# Constructor for single-run NMF results

Stores the W (genes x components) and H (components x cells) matrices
from a single HALS NMF run, plus convergence metadata and provenance
information.

## Usage

``` r
new_nmf_result(nmf_res, gene_ids, cell_ids, cell_indices, source_class, params)
```

## Arguments

- nmf_res:

  Named list. Output of
  [`rs_nmf_single_sc()`](https://gregorlueg.github.io/bixverse/reference/rs_nmf_single_sc.md)
  or
  [`rs_nmf_single_mc()`](https://gregorlueg.github.io/bixverse/reference/rs_nmf_single_mc.md).
  Must contain `w`, `h`, `final_loss`, `n_iter`, `converged`.

- gene_ids:

  Character vector. Gene identifiers for the rows of W.

- cell_ids:

  Character vector. Cell (or meta cell) identifiers for the columns of
  H.

- cell_indices:

  Integer vector. 0-indexed cell positions used in the run (Rust-side
  indices, kept for re-running / cross-referencing).

- source_class:

  String. One of `c("SingleCells", "MetaCells")`.

- params:

  List. The full set of parameters used for the run.

## Value

An object of class `NmfResult`.
