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

## Examples

``` r
# mean expression and expressing fraction per planted cell type
sc <- demo_single_cells()
res <- compute_expression_info_sc(
  sc,
  celltype_colname = "cell_grp",
  genes = get_gene_names(sc)[1:5]
)
head(res)
#>     cluster_id    gene avg_expr frac_expr
#>         <char>  <char>    <num>     <num>
#> 1: cell_type_1 gene_01 6.241650 0.9880240
#> 2: cell_type_1 gene_02 6.115141 0.9880240
#> 3: cell_type_1 gene_03 6.385783 1.0000000
#> 4: cell_type_1 gene_04 5.770888 0.9580838
#> 5: cell_type_1 gene_05 4.425770 0.8383234
#> 6: cell_type_2 gene_01 2.902391 0.6946108

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
