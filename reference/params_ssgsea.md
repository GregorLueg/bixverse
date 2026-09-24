# Wrapper function to generate ssGSEA parameters

Wrapper function to generate ssGSEA parameters

## Usage

``` r
params_ssgsea(alpha = 0.25, min_size = 5L, max_size = 500L, normalise = TRUE)
```

## Arguments

- alpha:

  Numeric. The exponent defining the weight of the tail in the random
  walk performed by ssGSEA. Defaults to `0.25`.

- min_size:

  Integer. Minimum number of genes per gene set. Defaults to `5L`.

- max_size:

  Integer. Maximum number of genes per gene set. Defaults to `500L`.

- normalise:

  Boolean. Shall the scores be normalised. Defaults to `TRUE`.

## Value

A named list with the following elements:

- alpha - Numeric. The exponent defining the weight of the tail in the
  random walk performed by ssGSEA. Defaults to `0.25`.

- min_size - Integer. Minimum number of genes per gene set. Defaults to
  `5L`.

- max_size - Integer. Maximum number of genes per gene set. Defaults to
  `500L`.

- normalise - Boolean. Shall the scores be normalised. Defaults to
  `TRUE`.
