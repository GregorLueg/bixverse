# Similarity network fusion

This function iteratively fuses the affinity matrices together.

## Usage

``` r
rs_snf(aff_mat_list, k, t, alpha)
```

## Arguments

- aff_mat_list:

  A list of numerical matrices. The affinity matrices to fuse together.

- k:

  Integer. Number of neighbours to consider.

- t:

  Integer. Number of iterations for the algorithm.

- alpha:

  Float. Normalisation parameter controlling the fusion strength.

## Value

The final affinity matrix after the fusion.
