# Melt a PAGA connectivity matrix into a one-row-per-edge table

Both abstracted graphs are stored symmetrically with a zero diagonal, so
the lower triangle is dropped here. Keeping it would draw every edge
twice, at twice the apparent width.

## Usage

``` r
.paga_edges(conn, threshold, keep)
```

## Arguments

- conn:

  Sparse matrix. The abstracted graph, named by cluster.

- threshold:

  Numeric. Edges below this connectivity are dropped.

- keep:

  Character vector. Clusters that survived the empty-cluster filter.
  Edges touching anything else are dropped, as they have no end point to
  attach to.

## Value

A data.table with `from`, `to` and `weight`.
