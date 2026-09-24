# Calculates PCA on Pearson residuals for single cell

**\[experimental\]** Runs the PCA on the residuals a fitted model
implies, rather than on the stored normalised layer. The residual
columns are dense by construction, since a zero count still has a
residual, so there is no sparse or streaming variant of this path.

Two settings are refused rather than ignored: the `PFlogPF` transform,
which belongs to the normalised layer, and variance normalisation, which
would flatten the very ranking the residuals produce.

## Usage

``` r
rs_sc_pca_residuals(
  f_path_gene,
  residual_fit,
  no_pcs,
  pca_params,
  cell_indices,
  gene_indices,
  seed,
  return_scaled,
  verbose
)
```

## Arguments

- f_path_gene:

  String. Path to the `counts_genes.bin` file.

- residual_fit:

  List. A fit from
  [`rs_sc_fit_residuals()`](https://gregorlueg.github.io/bixverse/reference/rs_sc_fit_residuals.md).

- no_pcs:

  Integer. Number of PCs to calculate.

- pca_params:

  Named list. Contains the parameters to use for this PCA run. `clr` and
  `normalise_variance` must both be `FALSE`.

- cell_indices:

  Integer vector. The cell indices to use. (0-indexed!) Must be the
  selection the fit was fitted on.

- gene_indices:

  Integer vector. The gene indices to use. (0-indexed!) Every one must
  be covered by the fit.

- seed:

  Integer. Random seed for the randomised SVD.

- return_scaled:

  Boolean. Shall the scaled data be returned.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following items

- scores - The samples projected on the PCA space.

- loadings - The loadings of the features for the PCA.

- singular_values - The singular values for the PCA.

- scaled - The scaled matrix if `return_scaled = TRUE`, otherwise
  `NULL`.
