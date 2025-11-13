# Get the contrastive PCA factors

Getter function for the feature factors of the contrastive PCA

## Usage

``` r
get_c_pca_factors(object)
```

## Arguments

- object:

  The underlying class, see [`bulk_coexp()`](bulk_coexp.md).

## Value

The sample scores of the contrastive PCA run. If not found, returns a
warning and NULL.
