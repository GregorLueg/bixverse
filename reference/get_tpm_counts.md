# Return the TPM-normalised counts

Getter function to extract the TPM-normalised counts from the
[`BulkDge()`](https://gregorlueg.github.io/bixverse/reference/BulkDge.md)
class.

## Usage

``` r
get_tpm_counts(object)
```

## Arguments

- object:

  `BulkDge` class.

## Value

Returns the TPM-normalised counts. (If found.)

## Examples

``` r
# TPM counts, which normalise_bulk_dge() only computes on request
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
  calc_tpm = TRUE,
  gene_lengths = gene_lengths,
  .verbose = FALSE
)
get_tpm_counts(object)[1:3, 1:3]
#>         sample_1 sample_10 sample_100
#> gene_1 726.28034 498.21473   506.6152
#> gene_2   0.00000   0.00000    77.9408
#> gene_3  44.01699  62.27684   759.9228
```
