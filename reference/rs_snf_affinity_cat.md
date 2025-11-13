# Calculate the SNF affinity matrix for categorical values

Calculate the SNF affinity matrix for categorical values

## Usage

``` r
rs_snf_affinity_cat(data, k, mu)
```

## Arguments

- data:

  Integer matrix. Needs to be oriented features x samples! The integers
  represent the factor values of the catagories.

- k:

  Integer. Number of neighbours to consider.

- mu:

  Float. Normalisation factor for the Gaussian kernel width.

## Value

The affinity matrix based on categorical values.
