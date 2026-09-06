# Return the outputs

Getter function to extract the outputs from the
[`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md)
or
[`BulkDge()`](https://gregorlueg.github.io/bixverse/reference/BulkDge.md).

## Usage

``` r
get_outputs(object, ...)
```

## Arguments

- object:

  The underlying object, either
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md)
  or
  [`BulkDge()`](https://gregorlueg.github.io/bixverse/reference/BulkDge.md).

- ...:

  Additional arguments to parse to the functions.

## Value

Returns the outputs stored in the class.

## Examples

``` r
# everything the QC and normalisation steps stashed on the object
syn <- synthetic_bulk_cor_matrix()
meta <- data.table::data.table(
  sample_id = colnames(syn$counts),
  case_control = rep(c("case", "control"), each = 50)
)
object <- BulkDge(raw_counts = syn$counts, meta_data = meta)
object <- qc_bulk_dge(object, group_col = "case_control", .verbose = FALSE)
names(get_outputs(object))
#> [1] "dge_list"            "sample_info"         "group_col"          
#> [4] "raw_counts_filtered"
```
