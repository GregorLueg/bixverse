# Prepare class for contrastive PCA

This function will prepare the `BulkCoExp` for subsequent usage of the
contrastive PCA functions. This is based on the work of Abid, et al.

## Usage

``` r
contrastive_pca_processing(
  object,
  background_matrix,
  scale = FALSE,
  .verbose = TRUE
)
```

## Arguments

- object:

  The underlying class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

- background_matrix:

  Numeric matrix. The background matrix you wish to remove. You should
  apply any data transformation to this matrix, too!

- scale:

  Boolean. Shall the data be scaled. Defaults to FALSE.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

`BulkCoExp` with the needed data for contrastive PCA in the properties
of the class.

## References

Abid, et al., Nature Communications, 2018

## Examples

``` r
# covariance matrices of the target and the background
cpca_data <- synthetic_c_pca_data()
target <- t(cpca_data$target)
background <- t(cpca_data$background)
meta <- data.table::data.table(sample_id = rownames(target))
obj <- BulkCoExp(target, meta)
obj <- preprocess_bulk_coexp(obj, .verbose = FALSE)
obj <- contrastive_pca_processing(obj, background, .verbose = FALSE)
dim(obj@processed_data$target_covar)
#> [1] 30 30
```
