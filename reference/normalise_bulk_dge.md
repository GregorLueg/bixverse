# Normalise the count data for DGE.

This function will apply the CPM + Voom normalisation and can
additionally calculate TPM and FPKM values for plotting purposes.

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
  Please refer to
  [`edgeR::normLibSizes()`](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html).

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
#> calcNormFactors has been renamed to normLibSizes
get_outputs(object)$normalised_counts[1:3, 1:3]
#>        sample_1 sample_10 sample_100
#> gene_1 9.529540  8.997354   8.994702
#> gene_2 4.848276  4.259485   6.473063
#> gene_3 6.473063  6.061875   9.555580
```
