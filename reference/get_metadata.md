# Return the metadata

Getter function to extract the metadata from the
[`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md)
or
[`BulkDge()`](https://gregorlueg.github.io/bixverse/reference/BulkDge.md).

## Usage

``` r
get_metadata(object, ...)
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

Returns the metadata stored in the class.

## Examples

``` r
# sample metadata back out of a BulkDge
syn <- synthetic_bulk_cor_matrix()
meta <- data.table::data.table(
  sample_id = colnames(syn$counts),
  case_control = rep(c("case", "control"), each = 50)
)
object <- BulkDge(raw_counts = syn$counts, meta_data = meta)
head(get_metadata(object))
#>    sample_id case_control
#>       <char>       <char>
#> 1:  sample_1         case
#> 2:  sample_2         case
#> 3:  sample_3         case
#> 4:  sample_4         case
#> 5:  sample_5         case
#> 6:  sample_6         case
```
