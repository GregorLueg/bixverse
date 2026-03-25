# Plot per-cell QC violin plots from a CellQc object

Plot per-cell QC violin plots from a CellQc object

## Usage

``` r
# S3 method for class 'CellQc'
plot(x, qc_df, ...)
```

## Arguments

- x:

  A `CellQc` object.

- qc_df:

  A data.table containing the cell-level data.

- ...:

  Ignored.

## Value

A named list of ggplot objects, one per metric.
