# Constructor for consensus NMF results

Stores the consensus W (genes x components) and H (components x cells)
matrices, the relative reconstruction errors and the clustering
diagnostics that produced them.

## Usage

``` r
new_consensus_nmf_result(
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
  [`rs_nmf_consensus_sc()`](https://gregorlueg.github.io/bixverse/reference/rs_nmf_consensus_sc.md)
  or
  [`rs_nmf_consensus_mc()`](https://gregorlueg.github.io/bixverse/reference/rs_nmf_consensus_mc.md).
  Must contain `w`, `h`, `rel_error`, `rel_run_errors`, `labels`,
  `local_density`, `kept`, `silhouette`, `stability`, `cluster_sizes`,
  `n_dropped`, `n_empty_clusters`.

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

An object of class `ConsensusNmfResult`.

## References

Kotliar et al., eLife, 2019
