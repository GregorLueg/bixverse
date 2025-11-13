# Calculates the Limma Voom DGE

This function will apply the Limma Voom DGE workflow. At a minimum you
will need to provide `contrast_column` that can be found in the
meta-data. If you do not provide a vector of contrasts that you wish to
test for, every permutation of groups represented in that column will be
tested against each other.

## Usage

``` r
calculate_dge_limma(
  object,
  contrast_column,
  contrast_list = NULL,
  filter_column = NULL,
  co_variates = NULL,
  quantile_norm = FALSE,
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
  tissue separately in the data.

- co_variates:

  Optional string vector. Any co-variates you wish to consider during
  the Limma Voom modelling.

- quantile_norm:

  Boolean. Shall the data also be quantile normalised. Defaults to
  `FALSE`.

- ...:

  Additional parameters to forward to
  [`limma::eBayes()`](https://rdrr.io/pkg/limma/man/ebayes.html) or
  [`limma::voom()`](https://rdrr.io/pkg/limma/man/voom.html).

- .verbose:

  Controls verbosity of the function.

## Value

Returns the class with additional data added to the outputs.
