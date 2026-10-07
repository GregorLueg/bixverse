# Helper function to generate default parameters for PCA

Helper function to generate default parameters for PCA

## Usage

``` r
params_pca_defaults()
```

## Value

A named list with the following elements:

- no_pcs - Integer. Number of PCs to consider. Defaults to `30L`.

- svd_solver - String. Which solver to use. `"randomised"` (default) is
  a randomised SVD, approximate in the trailing components.
  `"covariance"` builds the gene x gene cross-product and
  eigendecomposes it. `"exact"` is Lanczos on the sparse path and a full
  SVD on the dense one. See
  [`params_sc_pca()`](https://gregorlueg.github.io/bixverse/reference/params_sc_pca.md)
  for the trade-offs. One of `c("randomised", "covariance", "exact")`.
  Defaults to `"randomised"`.
