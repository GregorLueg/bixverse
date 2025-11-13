# Generates CoReMo-based gene modules

This function creates gene modules, based on the framework from
Srivastava et al., 2018. Briefly, it applies an RBF function to the
correlation matrix to reduce the impact of weak correlations, leverages
hierarchical clustering for clustering the genes. It optimises the R2
(cor^2) within the clusters to identify the optimal cut. Gene modules
with low R2 are being considered as the 'junk module'.

## Usage

``` r
cor_module_coremo_clustering(
  object,
  span = 0.1,
  coremo_params = params_coremo(),
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  The class, see [`bulk_coexp()`](bulk_coexp.md).

- span:

  Float. Defines the span parameter for the loess function to identify
  the inflection point. Defaults to `0.1`. Must be between 0.1 and 1.

- coremo_params:

  List. Parameters for the generation of the CoReMo modules, see
  [`params_coremo()`](params_coremo.md). Contains:

  - epsilon - Float. The epsilon parameter for the RBF. You can optimise
    that one with
    [`cor_module_check_epsilon()`](cor_module_check_epsilon.md).
    Defaults to `2`.

  - k_min - Integer. Minimum number of cuts. Defaults to `2L`.

  - k_max - Integer. Maximum number of cuts. Defaults to `150L`.

  - min_size - Optional integer. Minimum size of the clusters. If
    provided, smaller clusters will be merged by eigengene similarity.

  - junk_module_threshold - Float. Minimum R2 median value for a module
    to not be considered a junk module. Defaults to `0.05`.

  - rbf_func - String. Type of RBF function to apply. Defaults to
    `"gaussian"`.

  - cor_method - String. Type of correlation method to use for merging
    the smaller cluster. Defaults to `"spearman"`.

- seed:

  Integer. Random seed for reproducibility purposes.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

The class with added data to the properties for subsequent usage.

## References

Srivastava, et al., Nat. Commun., 2018
