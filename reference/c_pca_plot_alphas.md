# Plot various alphas for the contrastive PCA

This function will plot various alphas to highlight the most interesting
alpha parameters akin to the implementation of contrastive PCA in
Python.

## Usage

``` r
c_pca_plot_alphas(
  object,
  label_column = NULL,
  min_alpha = 0.1,
  max_alpha = 100,
  n_alphas = 10L,
  .verbose = TRUE
)
```

## Arguments

- object:

  The underlying class, see [`bulk_coexp()`](bulk_coexp.md). You need to
  apply [`contrastive_pca_processing()`](contrastive_pca_processing.md)
  to the function for this method to work. Checkmate will raise errors
  otherwise.

- label_column:

  An optional sample label column. Needs to exist in the meta_data of
  the `bulk_coexp` class.

- min_alpha:

  Minimum alpha to test.

- max_alpha:

  Maximum alpha to test.

- n_alphas:

  Number of alphas to test. The function will generate a series of
  alphas from log(min_alpha) to log(max_alpha) to test out.

- .verbose:

  Controls verbosity of function.

## Value

A ggplot showing the impact of various alpha parameters on the samples
in form of 2D plots.

## References

Abid, et al., Nature Communications, 2018
