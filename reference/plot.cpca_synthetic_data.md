# Plot the contrastive PCA example data

Two heatmaps side by side, the target and the background matrix, so you
can see the structure contrastive PCA is meant to pull apart.

## Usage

``` r
# S3 method for class 'cpca_synthetic_data'
plot(x, ...)
```

## Arguments

- x:

  `cpca_synthetic_data` class. Output from
  [`synthetic_c_pca_data()`](https://gregorlueg.github.io/bixverse/reference/synthetic_c_pca_data.md).

- ...:

  Additional params

## Value

A ggplot showing the two heatmaps from the target and background matrix.

## Examples

``` r
# target next to background, the structure cPCA pulls apart
cpca_data <- synthetic_c_pca_data()
plot(cpca_data)

```
