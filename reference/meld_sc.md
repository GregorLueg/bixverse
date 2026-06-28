# Run MELD signal smoothing for differential abundance estimation

This function implements MELD for estimating sample-associated density
on a cell manifold. The general idea is to smooth a binary sample
indicator matrix over the kNN graph via spectral filtering, yielding
per-cell likelihood estimates for each sample condition. For further
details on the method, please refer to Burkhardt, et al. This function
will take a `SingleCells` class and return the smoothed density
estimates per cell and condition.

## Usage

``` r
meld_sc(
  object,
  sample_id_col,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  meld_params = params_meld(),
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` class.

- sample_id_col:

  Character. The column in the obs table representing the sample
  identifier.

- embd_to_use:

  Character. The embedding to use for kNN graph construction. Please use
  the same here as you used to generate the neighbours. Defaults to
  `"pca"`.

- no_embd_to_use:

  Optional integer. If you only want to use a subset of the embedding
  dimensions.

- meld_params:

  A list, please see
  [`params_meld()`](https://gregorlueg.github.io/bixverse/reference/params_meld.md).
  The list has the following parameters:

  - beta - Numeric. Smoothing strength; larger values produce smoother
    densities.

  - offset - Numeric. Shift of the filter centre in the rescaled
    spectrum. Must be in `[0, 1]`.

  - order - Numeric. Filter falloff sharpness; larger values approach a
    square low-pass.

  - filter - Character. Filter family. One of `c("heat", "laplacian")`.

  - chebyshev_order - Integer. Number of Chebyshev coefficients. Must be
    \>= 2.

  - lap_type - Character. Type of Laplacian. One of
    `c("combinatorial", "normalised")`.

  - normalise_indicators - Logical. If `TRUE`, indicator columns are
    divided by their column sum before filtering.

  - landmark - Logical. Whether to use landmark approximation. Useful
    for large data sets.

  - n_landmarks - Integer. Number of landmarks to use.

  - knn - List of kNN parameters. See
    [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
    for available parameters and their defaults.

- seed:

  Integer. Seed for reproducibility.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

A list with:

- raw_scores - The raw MELD scores

- norm_scores - Negative values were clamped to 0 and the rows L1
  normalised. This yields probability-like values.

## References

Burkhardt, et al. Nat. Biotechnol., 2021.
