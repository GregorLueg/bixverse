# Calculate the contrastive PCA

**\[experimental\]** This function calculates the contrastive PCA given
a target covariance matrix and the background covariance matrix you wish
to subtract. The alpha parameter controls how much of the background
covariance you wish to remove. You have the options to return the
feature loadings and you can specify the number of cPCAs to return.

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

  The original values of the target matrix. Rows = samples, columns =
  features.

- alpha:

  How much of the background co-variance should be removed.

- n_pcs:

  How many contrastive PCs to return

- return_loadings:

  Shall the loadings be returned from the contrastive PCA

## Value

A list containing:

- factors - The factors of the contrastive PCA, i.e. `target_mat`
  multiplied by the loadings. Samples x `n_pcs`.

- loadings - The loadings (top eigenvectors) of the contrastive PCA.
  Features x `n_pcs`. Will be `NULL` if `return_loadings = FALSE`.
