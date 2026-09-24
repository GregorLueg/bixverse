# Residual-based PCA for meta cells

Residual-based PCA for meta cells

## Usage

``` r
.calculate_pca_residual_mc(
  object,
  no_pcs,
  pca_params,
  selected_hvg,
  seed,
  .verbose
)
```

## Arguments

- object:

  `MetaCells` class.

- no_pcs:

  Integer. Number of PCs.

- pca_params:

  List. See
  [`params_sc_pca()`](https://gregorlueg.github.io/bixverse/reference/params_sc_pca.md).

- selected_hvg:

  Integer. The 1-based genes to use.

- seed:

  Integer. Random seed.

- .verbose:

  Boolean or Integer. Controls verbosity.

## Value

The object with the PCA attached.
