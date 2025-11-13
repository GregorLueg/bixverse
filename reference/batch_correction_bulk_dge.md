# Run a linear batch correction

Runs a linear batch correction over the data regressing out batch
effects and adding `normalised_counts_corrected` to the object. Should
these counts be found by
[`calculate_dge_hedges()`](calculate_dge_hedges.md), they will be used
for calculations of effect sizes based on Hedge's G.

## Usage

``` r
batch_correction_bulk_dge(
  object,
  contrast_column,
  batch_col,
  scale_genes = FALSE,
  no_hvg_genes = 2500L
)
```

## Arguments

- object:

  The underlying class, see [`bulk_dge()`](bulk_dge.md).

- contrast_column:

  String. The contrast column in which the groupings are stored. Needs
  to be found in the meta_data within the properties.

- batch_col:

  String. The column in which the batch effect groups can be found.

- scale_genes:

  Boolean. Shall the log(cpm) counts be scaled prior the PCA
  calculation. Defaults to `FALSE`.

- no_hvg_genes:

  Integer. Number of highly variable genes to include. Defaults to 2500.

## Value

Returns the class with additional data added to the outputs.
