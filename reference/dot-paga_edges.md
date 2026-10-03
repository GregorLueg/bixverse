# Melt a symmetric graph matrix into a one-row-per-edge table

Used for the PAGA abstracted graphs and the miloR neighbourhood overlap.
Both are symmetric, so the lower triangle and the diagonal are dropped
here. Keeping it would draw every edge twice, at twice the apparent
width.

## Usage

``` r
.paga_edges(conn, threshold, keep)
```

## Arguments

- conn:

  Sparse matrix. The graph, named by node.

- threshold:

  Numeric. Edges below this connectivity are dropped.

- keep:

  Character vector. Nodes that survived upstream filtering. Edges
  touching anything else are dropped, as they have no end point to
  attach to.

## Value

A data.table with `from`, `to` and `weight`.
