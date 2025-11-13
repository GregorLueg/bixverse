# QC on the bulk dge data (DEPRECATED!)

This is a deprecated method and will raise an error. Please use
[`qc_bulk_dge()`](qc_bulk_dge.md) and
[`normalise_bulk_dge()`](normalise_bulk_dge.md) instead.

## Usage

``` r
preprocess_bulk_dge(
  object,
  group_col,
  norm_method = c("TMM", "TMMwsp", "RLE", "upperquartile", "none"),
  outlier_threshold = 2,
  min_prop = 0.2,
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

- outlier_threshold:

  Float. Number of standard deviations in terms of percentage genes
  detected you allow before removing a sample. Defaults to `2`.

- min_prop:

  Float. Minimum proportion of samples in which the gene has to be
  identified in.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

Throws an error, as it is deprecated.
