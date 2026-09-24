# Calculate the principal component regression on batch

Regresses each embedding dimension on the batch labels and weights the
per-dimension R-squared by the variance of that dimension. On its own
the number says how much of the embedding variance batch explains. For a
corrected embedding, the function also runs it on the uncorrected PCA
and reports the scIB comparison `(pre - post) / pre`: 1 means the batch
variance is gone, 0 means nothing changed, negative means it got worse.

## Usage

``` r
calculate_pcr_sc(object, batch_column, embd_to_use = "pca", .verbose = TRUE)
```

## Arguments

- object:

  `SingleCells` or `SingleCellsSubset` class.

- batch_column:

  String. The column with the batch information in the obs data of the
  class.

- embd_to_use:

  String. Which embedding to compute the PCR on. Defaults to `"pca"`.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

A `PcrScores` object with the following elements

- pcr - Variance-weighted R-squared of batch on `embd_to_use`.

- pcr_pca - The same on the uncorrected PCA.

- pcr_comparison - `(pcr_pca - pcr) / pcr_pca`. `NA` if
  `embd_to_use = "pca"`.

- var_explained - Variance per embedding dimension.

- r_squared - Batch R-squared per embedding dimension.

- embedding_used - Which embedding the PCR was computed on.

## References

Büttner, et al., Nat. Methods, 2019; Luecken, et al., Nat. Methods, 2022

## Examples

``` r
# share of PCA variance explained by batch
sc <- demo_single_cells(
  syn_data_params = params_sc_synthetic_data(
    n_cells = 600L, n_genes = 50L, n_batches = 3L
  )
)
calculate_pcr_sc(sc, batch_column = "batch_index")
#> Principal Component Regression (batch)
#>   Embedding: pca | Dimensions: 10
#>   PCR:             0.1414

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
