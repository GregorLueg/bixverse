# Calculates PCA on Pearson residuals for meta cells

**\[experimental\]** In-memory version of
[`rs_sc_pca_residuals()`](https://gregorlueg.github.io/bixverse/reference/rs_sc_pca_residuals.md).
As there, `clr` and `normalise_variance` must both be `FALSE`.

## Usage

``` r
rs_mc_pca_residuals(
  sparse_data,
  residual_fit,
  no_pcs,
  pca_params,
  gene_indices,
  seed,
  verbose
)
```

## Arguments

- sparse_data:

  A named list that needs to have `data`, `indptr`, `indices`, `nrow`,
  `ncol` and `cs_type`. Shape is (metacells, genes). Pass raw counts.

- residual_fit:

  List. A fit from
  [`rs_mc_fit_residuals()`](https://gregorlueg.github.io/bixverse/reference/rs_mc_fit_residuals.md).

- no_pcs:

  Integer. Number of PCs to calculate.

- pca_params:

  Named list. Contains the parameters to use for this PCA run.

- gene_indices:

  Integer vector. The gene indices to use. (0-indexed!)

- seed:

  Integer. Random seed for the randomised SVD.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following items

- scores - The samples projected on the PCA space.

- loadings - The loadings of the features for the PCA.

- singular_values - The singular values for the PCA.
