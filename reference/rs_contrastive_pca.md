# Calculate the contrastive PCA

This function calculate the contrastive PCA given a target covariance
matrix and the background covariance matrix you wish to subtract. The
alpha parameter controls how much of the background covariance you wish
to remove. You have the options to return the feature loadings and you
can specificy the number of cPCAs to return. WARNING! Incorrect use can
cause kernel crashes. Wrapper around the Rust functions with type checks
are provided in the package.

## Usage

``` r
rs_contrastive_pca(
  target_covar,
  background_covar,
  target_mat,
  alpha,
  n_pcs,
  return_loadings
)
```

## Arguments

- target_covar:

  The co-variance matrix of the target data set.

- background_covar:

  The co-variance matrix of the background data set.

- target_mat:

  The original values of the target matrix.

- alpha:

  How much of the background co-variance should be removed.

- n_pcs:

  How many contrastive PCs to return

- return_loadings:

  Shall the loadings be returned from the contrastive PCA

## Value

A list containing:

- factors - The factors of the contrastive PCA.

- loadings - The loadings of the contrastive PCA. Will be NULL if
  return_loadings is set to FALSE.
