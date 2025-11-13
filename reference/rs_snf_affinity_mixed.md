# Calculate the SNF affinity matrix for mixed values

Calculate the SNF affinity matrix for mixed values

## Usage

``` r
rs_snf_affinity_mixed(data, is_cat, k, mu)
```

## Arguments

- data:

  Numerical matrix. Needs to be oriented features x samples! This
  function will calculate the Gower distance under the hood for the
  affinity calculation.

- is_cat:

  Boolean vector. Which of the features are categorical. Needs to be of
  `nrow(data)`.

- k:

  Integer. Number of neighbours to consider.

- mu:

  Float. Normalisation factor for the Gaussian kernel width.

## Value

The affinity matrix based on mixed values.
