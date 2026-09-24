# Return the FPKM-normalised counts

Getter function to extract the FPKM-normalised counts from the
[`BulkDge()`](https://gregorlueg.github.io/bixverse/reference/BulkDge.md)
class.

## Usage

``` r
get_fpkm_counts(object)
```

## Arguments

- object:

  `BulkDge` class.

## Value

Returns the FPKM-normalised counts. (If found.)

## Examples

``` r
# FPKM counts, which normalise_bulk_dge() only computes on request
syn <- synthetic_bulk_cor_matrix()
meta <- data.table::data.table(
  sample_id = colnames(syn$counts),
  case_control = rep(c("case", "control"), each = 50)
)
gene_lengths <- stats::setNames(
  rep(2000, nrow(syn$counts)),
  rownames(syn$counts)
)
object <- BulkDge(raw_counts = syn$counts, meta_data = meta)
object <- qc_bulk_dge(object, group_col = "case_control", .verbose = FALSE)
object <- normalise_bulk_dge(
  object,
  group_col = "case_control",
  calc_fpkm = TRUE,
  gene_lengths = gene_lengths,
  .verbose = FALSE
)
get_fpkm_counts(object)[1:3, 1:3]
#>        sample_1 sample_10 sample_100
#> gene_1 363.1402 249.10737   253.3076
#> gene_2   0.0000   0.00000    38.9704
#> gene_3  22.0085  31.13842   379.9614
```
