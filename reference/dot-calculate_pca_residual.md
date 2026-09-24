# Residual-based PCA

Shared body of the residual branch in
[`calculate_pca_sc()`](https://gregorlueg.github.io/bixverse/reference/calculate_pca_sc.md).
Enforces the settings the Rust side refuses, since those errors name
internals rather than the argument to change.

## Usage

``` r
.calculate_pca_residual(
  object,
  no_pcs,
  pca_params,
  selected_hvg,
  sparse_svd,
  seed,
  .verbose
)
```

## Arguments

- object:

  `SingleCells` or `SingleCellsSubset` class.

- no_pcs:

  Integer. Number of PCs.

- pca_params:

  List. See
  [`params_sc_pca()`](https://gregorlueg.github.io/bixverse/reference/params_sc_pca.md).

- selected_hvg:

  Integer. The 0-based genes to use.

- sparse_svd:

  Boolean. Must be `FALSE`.

- seed:

  Integer. Random seed.

- .verbose:

  Boolean or Integer. Controls verbosity.

## Value

The object with the PCA attached.
