# Apply contrastive PCA.

Applies the contrastive PCA algorithm given a specified alpha and a
number of contrastive principal components to extract.

## Usage

``` r
contrastive_pca(object, alpha, no_pcs)
```

## Arguments

- object:

  The underlying class, see [`bulk_coexp()`](bulk_coexp.md).

- alpha:

  Alpha parameter to use.

- no_pcs:

  Number of contrastive PCs to generate.

## Value

`bulk_coexp` with additional data in the slots

## References

Abid, et al., Nature Communications, 2018
