# Carry the input dimnames over to a missing-value PCA result

Carry the input dimnames over to a missing-value PCA result

## Usage

``` r
.name_missing_pca_res(res, x)
```

## Arguments

- res:

  List. Output of
  [`rs_ppca()`](https://gregorlueg.github.io/bixverse/reference/rs_ppca.md)
  or
  [`rs_bpca()`](https://gregorlueg.github.io/bixverse/reference/rs_bpca.md).

- x:

  Numeric matrix. The input matrix.

## Value

`res` with row names on `scores`, row names on `loadings` and the
dimnames of `x` on `completed`. Components are named `PC1`, `PC2`, etc.
