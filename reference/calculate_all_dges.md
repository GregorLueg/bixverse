# Calculate all possible DGE variants (DEPRECATED!)

This is a deprecated method and will raise an error. Please use
[`calculate_dge_limma()`](calculate_dge_limma.md) and
[`calculate_dge_hedges()`](calculate_dge_hedges.md).

## Usage

``` r
calculate_all_dges(
  object,
  contrast_column,
  contrast_list = NULL,
  filter_column = NULL,
  co_variates = NULL,
  small_sample_correction = NULL,
  ...,
  .verbose = TRUE
)
```

## Arguments

- object:

  The underlying class, see [`bulk_dge()`](bulk_dge.md).

- contrast_column:

  String. The contrast column in which the groupings are stored. Needs
  to be found in the meta_data within the properties.

- contrast_list:

  Optional string vector. A vectors that contains the contrast formatted
  as `"contrast1-contrast2"`. Default `NULL` will create all possible
  contrast automatically.

- filter_column:

  Optional String. If there is a column you wish to use as sub
  groupings, this can be provided here. An example could be different
  sampled tissues and you wish to run the DGE analyses within each
  tissue separately.

- co_variates:

  Optional string vector. Any co-variates you wish to consider during
  the Limma Voom modelling.

- small_sample_correction:

  Can be NULL (automatic determination if a small sample size correction
  should be applied) or a Boolean.

- ...:

  Additional parameters to forward to
  [`limma::eBayes()`](https://rdrr.io/pkg/limma/man/ebayes.html) or
  [`limma::voom()`](https://rdrr.io/pkg/limma/man/voom.html).

- .verbose:

  Controls verbosity of the function.

## Value

Throws an error, as it is deprecated.
