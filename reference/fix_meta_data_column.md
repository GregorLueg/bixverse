# Helper to fix meta-data columns to be R conform

This function will update the specified columns in the metadata of an
[`BulkDge()`](https://gregorlueg.github.io/bixverse/reference/BulkDge.md)
or
[`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md)
to be conform with R standard naming conventions. This is useful to do
before running DGE methods as they expect standardised names.

## Usage

``` r
fix_meta_data_column(object, col_names, ...)
```

## Arguments

- object:

  The underlying object, either
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md)
  or
  [`BulkDge()`](https://gregorlueg.github.io/bixverse/reference/BulkDge.md).

- col_names:

  Character vector. The columns to fix.

- ...:

  Additional arguments to parse to the functions.

## Value

Returns the object with the respective metadata columns updated.

## Examples

``` r
# make a contrast column safe for model.matrix()
syn <- synthetic_bulk_cor_matrix()
meta <- data.table::data.table(
  sample_id = colnames(syn$counts),
  case_control = rep(c("case 1", "control 1"), each = 50)
)
object <- BulkDge(raw_counts = syn$counts, meta_data = meta)
object <- fix_meta_data_column(object, "case_control")
unique(get_metadata(object)$case_control)
#> [1] case_1    control_1
#> Levels: case_1 control_1
```
