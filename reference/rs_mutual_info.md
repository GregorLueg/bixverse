# Calculates the mutual information matrix

**\[experimental\]** Calculates the pairwise mutual information across
all columns in the data.

## Usage

``` r
rs_mutual_info(x, n_bins, strategy, normalise)
```

## Arguments

- x:

  R matrix with doubles for which to calculate the mutual information

- n_bins:

  Optional integer. Number of bins to use. If `NULL` is provided the
  function will default to `sqrt(nrow(x))`.

- strategy:

  String. Binning strategy. One of `c("equal_width", "equal_freq")`.
  Unknown strings default to `"equal_width"`.

- normalise:

  Boolean. Shall the normalised mutual information be calculated via
  joint entropy.

## Value

The symmetric mutual information matrix. The diagonal holds the column
entropies, or `0` if `normalise = TRUE`.
