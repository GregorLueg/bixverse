# Get the contrastive PCA loadings

Getter function for the feature loadings of the contrastive PCA

## Usage

``` r
get_c_pca_loadings(object)
```

## Arguments

- object:

  The underlying class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

## Value

The feature loadings of the contrastive PCA run. If not found, returns a
warning and NULL.

## Examples

``` r
# feature loadings of the contrastive components
cpca_data <- synthetic_c_pca_data()
target <- t(cpca_data$target)
background <- t(cpca_data$background)
meta <- data.table::data.table(sample_id = rownames(target))
obj <- BulkCoExp(target, meta)
obj <- preprocess_bulk_coexp(obj, .verbose = FALSE)
obj <- contrastive_pca_processing(obj, background, .verbose = FALSE)
obj <- contrastive_pca(obj, alpha = 2.5, no_pcs = 5L)
head(get_c_pca_loadings(obj)[, 1:3])
#>                   cPC_1        cPC_2         cPC_3
#> feature_1  -0.020646762  0.030973018 -0.0133776245
#> feature_9   0.034373886  0.005184735  0.0070808801
#> feature_10  0.018730361  0.001767020 -0.0008550033
#> feature_5  -0.008601644  0.010760683 -0.0101120801
#> feature_8  -0.001500864  0.009470081 -0.0230014453
#> feature_2   0.003709886 -0.013213907  0.0043487894
```
