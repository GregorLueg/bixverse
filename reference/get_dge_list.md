# Return the DGEList

Builds an edgeR `DGEList` from the filtered counts, library sizes and
normalisation factors stored in the
[`BulkDge()`](https://gregorlueg.github.io/bixverse/reference/BulkDge.md)
class. bixverse does not need edgeR itself, so this requires edgeR to be
installed. Normalisation factors are one until
[`normalise_bulk_dge()`](https://gregorlueg.github.io/bixverse/reference/normalise_bulk_dge.md)
has run.

## Usage

``` r
get_dge_list(object)
```

## Arguments

- object:

  `BulkDge` class.

## Value

An edgeR `DGEList`, or `NULL` with a warning if
[`qc_bulk_dge()`](https://gregorlueg.github.io/bixverse/reference/qc_bulk_dge.md)
has not been run.

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
if (requireNamespace("edgeR", quietly = TRUE)) {
  dim(get_dge_list(object))
}
#> [1] 997  98
```
