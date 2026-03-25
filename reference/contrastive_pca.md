# Apply contrastive PCA.

Applies the contrastive PCA algorithm given a specified alpha and a
number of contrastive principal components to extract.

## Usage

``` r
contrastive_pca(object, alpha, no_pcs)
```

## Arguments

- object:

  The underlying class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

- alpha:

  Alpha parameter to use.

- no_pcs:

  Number of contrastive PCs to generate.

## Value

`BulkCoExp` with additional data in the slots

## References

Abid, et al., Nature Communications, 2018
