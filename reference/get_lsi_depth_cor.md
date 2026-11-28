# Correlation of each LSI component with the library size

Correlation of each LSI component with the library size

## Usage

``` r
get_lsi_depth_cor(object)
```

## Arguments

- object:

  `SingleCellsAtac` or `SingleCellsMultiModal` class, after
  [`calculate_lsi_sc()`](https://gregorlueg.github.io/bixverse/reference/calculate_lsi_sc.md).

## Value

A data.table with `component` (1-based), `depth_cor` and `dropped`.
