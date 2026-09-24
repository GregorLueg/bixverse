# Normalise the count data for DGE.

Calculates the normalisation factors and applies voom on the filtered
counts from
[`qc_bulk_dge()`](https://gregorlueg.github.io/bixverse/reference/qc_bulk_dge.md),
both in Rust via the `edge-rs` crate. Can additionally calculate TPM and
FPKM values for plotting purposes.

## Usage

``` r
normalise_bulk_dge(
  object,
  group_col,
  norm_method = c("TMM", "TMMwsp", "RLE", "upperquartile", "none"),
  calc_tpm = FALSE,
  calc_fpkm = FALSE,
  gene_lengths = NULL,
  .verbose = TRUE
)
```

## Arguments

- object:

  The underlying class, see
  [`BulkDge()`](https://gregorlueg.github.io/bixverse/reference/BulkDge.md).

- group_col:

  String. The column in the metadata that will contain the contrast
  groups. Needs to be part of the metadata stored in the class.

- norm_method:

  String. One of `c("TMM", "TMMwsp", "RLE", "upperquartile", "none")`.
  Please refer to edgeR's `calcNormFactors()`.

- calc_tpm:

  Boolean. Output TPM calculation (default = FALSE).

- calc_fpkm:

  Boolean. Output FPKM calculation (default = FALSE).

- gene_lengths:

  Optional named numeric. If you want to calculate TPM or FPKM you need
  to provide this one. The names need to be the same identifier as used
  in the counts.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

Returns the class with the `processed_data` data slot populated and
applied parameters added to the `params` slot.

## Examples

``` r
# TMM library size normalisation followed by voom
syn <- synthetic_bulk_cor_matrix()
meta <- data.table::data.table(
  sample_id = colnames(syn$counts),
  case_control = rep(c("case", "control"), each = 50)
)
object <- BulkDge(raw_counts = syn$counts, meta_data = meta)
object <- qc_bulk_dge(object, group_col = "case_control", .verbose = FALSE)
object <- normalise_bulk_dge(
  object,
  group_col = "case_control",
  .verbose = FALSE
)
get_outputs(object)$normalised_counts[1:3, 1:3]
#>        sample_1 sample_10 sample_100
#> gene_1 9.438292  8.938609   8.962916
#> gene_2 3.372203  3.323899   6.404920
#> gene_3 5.694131  6.131254   9.538776
```
