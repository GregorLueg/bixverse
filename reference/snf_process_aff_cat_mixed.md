# Helper function to process categorical or mixed data for SNF

Helper function to process categorical or mixed data for SNF

## Usage

``` r
snf_process_aff_cat_mixed(data, k, mu)
```

## Arguments

- data:

  data.table of structure samples x features. The function will assume
  that the first column represents the sample identifiers.

- k:

  Integer. Number of neighbours to consider.

- mu:

  Float. Normalisation factor for the Gaussian kernel width.

## Value

The affinity matrix based on categorical and/or mixed values.
