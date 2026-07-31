# Constructor for the per-cell ScType results

Constructor for the per-cell ScType results

## Usage

``` r
new_sc_type_cell_results(res)
```

## Arguments

- res:

  List. The raw output of
  [`rs_sc_type_assign_cells()`](https://gregorlueg.github.io/bixverse/reference/rs_sc_type_assign_cells.md).

## Value

An `ScTypeCellResults` object, a list with

- cell_types - Character vector. The scored cell types.

- assignments - Character vector. Per-cell call, `NA` for Unknown.

- scores - Numeric vector. Winning score per cell.

- margins - Numeric vector. Best minus second best score per cell. A
  small margin flags an ambiguous call.

- agreement - Numeric vector. Fraction of sNN neighbours sharing the
  call. `NULL` if no graph was used.

- hybrid - Character vector. The hybrid call, `NULL` if no cluster
  column was provided.

- composition - A `data.table` with the per-cluster composition, or
  `NULL`.

- counts - Cluster by cell type count matrix, or `NULL`.
