# Calculate PCA on the expression.

Calculates the principal component on top of the filtered count matrix
and adds the information of the first 10 principal components to the
outputs.

## Usage

``` r
calculate_pca_bulk_dge(
  object,
  scale_genes = FALSE,
  pcs = 10L,
  no_hvg_genes = 2500L
)
```

## Arguments

- object:

  The underlying class, see [`bulk_dge()`](bulk_dge.md).

- scale_genes:

  Boolean. Shall the log(cpm) counts be scaled prior the PCA
  calculation. Defaults to `FALSE`.

- pcs:

  Integer. Number of PCs to return and add to the outputs slot.

- no_hvg_genes:

  Integer. Number of highly variable genes to include. Defaults to 2500.

## Value

Returns the class with additional data added to the outputs.
