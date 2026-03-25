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
