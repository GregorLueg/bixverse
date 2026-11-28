# Wrapper function for the ATAC LSI parameters

Parameters of
[`calculate_lsi_sc()`](https://gregorlueg.github.io/bixverse/reference/calculate_lsi_sc.md),
Signac's `RunSVD()` on the TF-IDF.

## Usage

``` r
params_atac_lsi(
  n_comp = 50L,
  drop_first = TRUE,
  solver = c("auto", "exact", "randomised", "covariance"),
  scale_max = NULL
)
```

## Arguments

- n_comp:

  Integer. Components to compute. Defaults to `50L`.

- drop_first:

  Boolean. Leave component 1 out of the embedding. It usually tracks
  sequencing depth;
  [`get_lsi_depth_cor()`](https://gregorlueg.github.io/bixverse/reference/get_lsi_depth_cor.md)
  shows whether it does. Defaults to `TRUE`.

- solver:

  String. SVD solver, see
  [`params_sc_pca()`](https://gregorlueg.github.io/bixverse/reference/params_sc_pca.md)
  for the trade-offs. One of
  `c("auto", "exact", "randomised", "covariance")`. Defaults to
  `"auto"`.

- scale_max:

  Numeric or `NULL`. Clip the z-scored components to `+-scale_max`.
  Defaults to `NULL`.

## Value

A named list with the following elements:

- n_comp - Integer. Components to compute. Defaults to `50L`.

- drop_first - Boolean. Leave component 1 out of the embedding. It
  usually tracks sequencing depth;
  [`get_lsi_depth_cor()`](https://gregorlueg.github.io/bixverse/reference/get_lsi_depth_cor.md)
  shows whether it does. Defaults to `TRUE`.

- solver - String. SVD solver, see
  [`params_sc_pca()`](https://gregorlueg.github.io/bixverse/reference/params_sc_pca.md)
  for the trade-offs. One of
  `c("auto", "exact", "randomised", "covariance")`. Defaults to
  `"auto"`.

- scale_max - Numeric or `NULL`. Clip the z-scored components to
  `+-scale_max`. Defaults to `NULL`.
