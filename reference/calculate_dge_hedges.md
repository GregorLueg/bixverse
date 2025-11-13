# Calculates the Hedge's G effect size

This function will calculate the Hedge's G effect size on the normalised
counts. Should batch-corrected counts be found, these will be used. At a
minimum you will need to provide `contrast_column` that can be found in
the meta-data. If you do not provide a vector of contrasts that you wish
to test for, every permutation of groups represented in that column will
be tested against each other.

## Usage

``` r
calculate_dge_hedges(
  object,
  contrast_column,
  contrast_list = NULL,
  filter_column = NULL,
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

- .verbose:

  Controls verbosity of the function.

## Value

Returns the class with additional data added to the outputs.
