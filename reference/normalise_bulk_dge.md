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

  The underlying class, see [`bulk_dge()`](bulk_dge.md).

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
