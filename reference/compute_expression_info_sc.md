# Compute per-cluster mean expression and expressing fraction for a gene set

Thin R wrapper around the Rust `compute_cluster_expression_stats`
routine. Streams gene chunks from the on-disk store and aggregates
expression across user-supplied cell clusters. Cells outside any cluster
are ignored.

If `condition_colname` and `condition_oi` are supplied, only cells from
that condition contribute to the aggregation.

## Usage

``` r
compute_expression_info_sc(
  object,
  celltype_colname,
  genes,
  condition_colname = NULL,
  condition_oi = NULL
)
```

## Arguments

- object:

  A `SingleCells` object.

- celltype_colname:

  Name of the cluster column in `obs`.

- genes:

  Character vector of gene IDs to aggregate over.

- condition_colname:

  Optional. Name of a condition column in `obs`.

- condition_oi:

  Optional. Value of `condition_colname` to subset to.

## Value

A long `data.table` with columns `cluster_id`, `gene`, `avg_expr`,
`frac_expr`.
