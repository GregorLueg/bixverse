# Return the DGEList

Getter function to extract the DGEList from the
[`BulkDge()`](https://gregorlueg.github.io/bixverse/reference/BulkDge.md)
class.

## Usage

``` r
get_dge_list(object)
```

## Arguments

- object:

  `BulkDge` class.

## Value

Returns the DGEList stored in the class.

## Examples

``` r
# the edgeR DGEList built during QC
syn <- synthetic_bulk_cor_matrix()
meta <- data.table::data.table(
  sample_id = colnames(syn$counts),
  case_control = rep(c("case", "control"), each = 50)
)
object <- BulkDge(raw_counts = syn$counts, meta_data = meta)
object <- qc_bulk_dge(object, group_col = "case_control", .verbose = FALSE)
dim(get_dge_list(object))
#> [1] 997  98
```
