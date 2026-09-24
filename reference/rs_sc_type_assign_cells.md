# Per-cell ScType assignment with optional kNN smoothing

**\[experimental\]** Assigns cell types per cell instead of per cluster.
If a graph is provided, the score matrix is smoothed via label spreading
over that graph before the per-cell argmax. If cluster labels are
provided, the per-cluster composition and the hybrid assignment (pure
clusters keep the cluster-level call, mixed clusters fall back to the
per-cell calls) are returned on top.

## Usage

``` r
rs_sc_type_assign_cells(sc_type_res, from, to, weights, cluster_labels, params)
```

## Arguments

- sc_type_res:

  List. The ScType results, see
  [`rs_sc_type()`](https://gregorlueg.github.io/bixverse/reference/rs_sc_type.md).

- from, to:

  Optional integer vectors. 1-indexed(!) edges of the sNN graph, each
  edge listed once; the graph is symmetrised on the way in. If either is
  `NULL`, no smoothing is applied.

- weights:

  Optional numeric vector. Edge weights, same length as `from`.
  Reciprocal weights are averaged. If `NULL`, the graph is unweighted.

- cluster_labels:

  Optional integer vector. 0-indexed(!) cluster assignment, of length of
  the scored cells.

- params:

  List. The output of
  [`params_sctype_cells()`](https://gregorlueg.github.io/bixverse/reference/params_sctype_cells.md).

## Value

A list with

- cell_types - String vector. The cell types.

- assignments - Integer vector. 1-based index into `cell_types` per
  cell, `0L` denoting Unknown.

- scores - Numeric vector. Winning score per cell.

- margins - Numeric vector. Best minus second best score per cell.

- agreement - Numeric vector. Fraction of graph neighbours sharing the
  call. `NULL` if no graph was provided.

- hybrid_assignments - Integer vector, as `assignments`. `NULL` if no
  `cluster_labels` were provided.

- composition - List with the per-cluster composition (`cluster_id`,
  `n_cells`, `n_unknown`, `dominant`, `second`, `purity`,
  `second_fraction`, `entropy`, `cluster_mixed` and the clusters x cell
  types count matrix `counts`). `dominant` and `second` use the same
  1-based encoding as `assignments`. `NULL` if no `cluster_labels` were
  provided.

## References

Ianevski et al., Nat Comm, 2022. Zhou et al., NIPS, 2004.
