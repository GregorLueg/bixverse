# Wrapper function to generate CoReMo parameters

Wrapper function to generate CoReMo parameters

## Usage

``` r
params_coremo(
  epsilon = 2,
  k_min = 2L,
  k_max = 150L,
  min_size = NULL,
  junk_module_threshold = 0.05,
  rbf_func = c("gaussian", "inverse_quadratic", "bump"),
  cor_method = c("spearman", "pearson")
)
```

## Arguments

- epsilon:

  Numeric. Epsilon parameter for the chosen RBF function, see
  `rbf_func`. The higher, the more aggressively low correlations will be
  shrunk. Defaults to `2.0`.

- k_min:

  Integer. Minimum and maximum number of cuts to use for the
  hierarchical clustering. Defaults to `2L`.

- k_max:

  Integer. Minimum and maximum number of cuts to use for the
  hierarchical clustering. Defaults to `150L`.

- min_size:

  Integer or `NULL`. Minimum size of the clusters. Smaller clusters will
  be combined together. Defaults to `NULL`.

- junk_module_threshold:

  Numeric. Threshold for the minimum correlation to be observed in a
  module. Defaults to `0.05`.

- rbf_func:

  String. Type of RBF you wish to apply to down-weigh weak correlations.
  One of `c("gaussian", "inverse_quadratic", "bump")`. Defaults to
  `"gaussian"`.

- cor_method:

  String. The type of correlation to use. One of
  `c("spearman", "pearson")`. Defaults to `"spearman"`.

## Value

A named list with the following elements:

- epsilon - Numeric. Epsilon parameter for the chosen RBF function, see
  `rbf_func`. The higher, the more aggressively low correlations will be
  shrunk. Defaults to `2.0`.

- k_min - Integer. Minimum and maximum number of cuts to use for the
  hierarchical clustering. Defaults to `2L`.

- k_max - Integer. Minimum and maximum number of cuts to use for the
  hierarchical clustering. Defaults to `150L`.

- min_size - Integer or `NULL`. Minimum size of the clusters. Smaller
  clusters will be combined together. Defaults to `NULL`.

- junk_module_threshold - Numeric. Threshold for the minimum correlation
  to be observed in a module. Defaults to `0.05`.

- rbf_func - String. Type of RBF you wish to apply to down-weigh weak
  correlations. One of `c("gaussian", "inverse_quadratic", "bump")`.
  Defaults to `"gaussian"`.

- cor_method - String. The type of correlation to use. One of
  `c("spearman", "pearson")`. Defaults to `"spearman"`.
