# Wrapper function to generate GSEA parameters

Wrapper function to generate GSEA parameters

## Usage

``` r
params_gsea(
  min_size = 5L,
  max_size = 500L,
  gsea_param = 1,
  sample_size = 101L,
  eps = 1e-50
)
```

## Arguments

- min_size:

  Integer. Minimum number of genes per gene set. Defaults to `5L`.

- max_size:

  Integer. Maximum number of genes per gene set. Defaults to `500L`.

- gsea_param:

  Numeric. GSEA parameter. Defaults to `1.0`.

- sample_size:

  Integer. Number of samples to iterate through for the multi-level
  implementation of fgsea. Defaults to `101L`.

- eps:

  Numeric. Boundary for calculating the p-value. Used for the multi-
  level implementation of fgsea. Defaults to `1e-50`.

## Value

A named list with the following elements:

- min_size - Integer. Minimum number of genes per gene set. Defaults to
  `5L`.

- max_size - Integer. Maximum number of genes per gene set. Defaults to
  `500L`.

- gsea_param - Numeric. GSEA parameter. Defaults to `1.0`.

- sample_size - Integer. Number of samples to iterate through for the
  multi-level implementation of fgsea. Defaults to `101L`.

- eps - Numeric. Boundary for calculating the p-value. Used for the
  multi- level implementation of fgsea. Defaults to `1e-50`.
