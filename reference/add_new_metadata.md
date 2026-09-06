# Replace the meta data

This function will replace the meta data within the given object

## Usage

``` r
add_new_metadata(object, new_metadata, ...)
```

## Arguments

- object:

  The class

- new_metadata:

  data.table. The new meta data you wish to add.

- ...:

  Additional arguments to parse to the functions.

## Value

The object with updated metadata.

## Examples

``` r
# swap in a metadata table that carries an extra batch column
set.seed(42)
counts <- matrix(rpois(60, 20), nrow = 10, ncol = 6)
rownames(counts) <- sprintf("gene_%i", 1:10)
colnames(counts) <- sprintf("sample_%i", 1:6)
meta <- data.table::data.table(
  sample_id = colnames(counts),
  case_control = rep(c("case", "control"), each = 3)
)
object <- BulkDge(raw_counts = counts, meta_data = meta)
new_meta <- data.table::copy(meta)[, batch := rep(c("b1", "b2"), 3)]
object <- add_new_metadata(object, new_metadata = new_meta)
head(S7::prop(object, "meta_data"))
#>    sample_id case_control  batch
#>       <char>       <char> <char>
#> 1:  sample_1         case     b1
#> 2:  sample_2         case     b2
#> 3:  sample_3         case     b1
#> 4:  sample_4      control     b2
#> 5:  sample_5      control     b1
#> 6:  sample_6      control     b2
```
