# Assign ScType cell types per cell, smoothing over a Rust neighbour graph

**\[experimental\]** As
[`rs_sc_type_assign_cells()`](https://gregorlueg.github.io/bixverse/reference/rs_sc_type_assign_cells.md),
but the graph is the Rust neighbour graph behind an external pointer
rather than an edge list.

## Usage

``` r
rs_sc_type_assign_cells_graph(sc_type_res, graph, cluster_labels, params)
```

## Arguments

- sc_type_res:

  List. The ScType results, see
  [`rs_sc_type()`](https://gregorlueg.github.io/bixverse/reference/rs_sc_type.md).

- graph:

  External pointer to the Rust neighbour graph, with one node per scored
  cell.

- cluster_labels:

  Optional integer vector. 0-indexed(!) cluster assignment, of length of
  the scored cells.

- params:

  List. The output of
  [`params_sctype_cells()`](https://gregorlueg.github.io/bixverse/reference/params_sctype_cells.md).

## Value

See
[`rs_sc_type_assign_cells()`](https://gregorlueg.github.io/bixverse/reference/rs_sc_type_assign_cells.md).

## References

Ianevski et al., Nat Comm, 2022. Zhou et al., NIPS, 2004.
