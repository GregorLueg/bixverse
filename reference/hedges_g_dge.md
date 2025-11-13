# Calculate the effect size

Calculate the effect size

## Usage

``` r
hedges_g_dge(
  meta_data,
  main_contrast,
  normalised_counts,
  contrast_list = NULL,
  .verbose = TRUE
)
```

## Arguments

- meta_data:

  data.table. The meta information about the experiment in which the
  contrast info can be found.

- main_contrast:

  String. Which column contains the main groups you want to calculate
  the Hedge's G effect for. Every permutation of the groups will be
  tested if `contrast_list` is `NULL`.

- normalised_counts:

  Numeric Matrix. The normalized count matrix.

- contrast_list:

  String vector or NULL. Optional string vector of contrast formatted as
  `"contrast1-contrast2"`. Default NULL will create all contrasts
  automatically.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

A data.table with the effect sizes and standard errors based on the
Hedge's G effect size for the groups.
