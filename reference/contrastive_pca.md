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

## Examples

``` r
# five contrastive PCs at alpha 2.5
cpca_data <- synthetic_c_pca_data()
target <- t(cpca_data$target)
background <- t(cpca_data$background)
meta <- data.table::data.table(sample_id = rownames(target))
obj <- BulkCoExp(target, meta)
obj <- preprocess_bulk_coexp(obj, .verbose = FALSE)
obj <- contrastive_pca_processing(obj, background, .verbose = FALSE)
obj <- contrastive_pca(obj, alpha = 2.5, no_pcs = 5L)
dim(get_c_pca_factors(obj))
#> [1] 400   5
```
