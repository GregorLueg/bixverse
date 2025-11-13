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

  Float. Epsilon parameter for the chosen RBF function, see `rbf_func`.
  The higher, the more aggressively low correlations will be shrunk.

- k_min, k_max:

  Integer. Minimum and maximum number of cuts to use for the
  hierarchical clustering.

- min_size:

  Optional integer. Minimum size of the clusters. Smaller clusters will
  be combined together.

- junk_module_threshold:

  Float. Threshold for the minimum correlation to be observed in a
  module. Defaults to `0.05`.

- rbf_func:

  String. Type of RBF you wish to apply to down-weigh weak correlations.
  Defaults to `"gaussian"`.

- cor_method:

  String. The type of correlation to use. Defaults to `"spearman"`.

## Value

List with parameters for usage in subsequent function.
