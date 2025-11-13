# Get the contrastive PCA loadings

Getter function for the feature loadings of the contrastive PCA

## Usage

``` r
get_c_pca_loadings(object)
```

## Arguments

- object:

  The underlying class, see [`bulk_coexp()`](bulk_coexp.md).

## Value

The feature loadings of the contrastive PCA run. If not found, returns a
warning and NULL.
