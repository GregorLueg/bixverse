# Wrapper function to generate ssGSEA parameters

Wrapper function to generate ssGSEA parameters

## Usage

``` r
params_ssgsea(alpha = 0.25, min_size = 5L, max_size = 500L, normalise = TRUE)
```

## Arguments

- alpha:

  Float. The exponent defining the weight of the tail in the random walk
  performed by ssGSEA.

- min_size:

  Integer. Minimum number of genes per gene set.

- max_size:

  Integer. Maximum number of genes per gene set.

- normalise:

  Boolean. Shall the scores be normalised.

## Value

List with parameters for usage in subsequent function.
