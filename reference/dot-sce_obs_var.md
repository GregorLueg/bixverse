# Extract the obs and var tables from a SingleCellExperiment

`colData` becomes obs and `rowData` becomes var, with the identifiers
put first because the first column of each is what becomes `cell_id` and
`gene_id` downstream.

Non-atomic columns are dropped. `colData` is a `DataFrame` and can carry
nested `DataFrame` or list columns, which have nowhere to go in the
DuckDB.

## Usage

``` r
.sce_obs_var(sce)
```

## Arguments

- sce:

  `SingleCellExperiment` class.

## Value

A list with `obs` and `var` as data.tables.
