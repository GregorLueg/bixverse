# Run a linear batch correction

Runs a linear batch correction over the data regressing out batch
effects and adding `normalised_counts_corrected` to the object. Should
these counts be found by
[`calculate_dge_hedges()`](https://gregorlueg.github.io/bixverse/reference/calculate_dge_hedges.md),
they will be used for calculations of effect sizes based on Hedge's G.

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

  The underlying class, see
  [`BulkDge()`](https://gregorlueg.github.io/bixverse/reference/BulkDge.md).

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

## Examples

``` r
# regress out a batch column and keep the corrected counts
syn <- synthetic_bulk_cor_matrix()
meta <- data.table::data.table(
  sample_id = colnames(syn$counts),
  case_control = rep(c("case", "control"), each = 50),
  batch = rep(c("b1", "b2"), times = 50)
)
object <- BulkDge(raw_counts = syn$counts, meta_data = meta)
object <- qc_bulk_dge(object, group_col = "case_control", .verbose = FALSE)
object <- normalise_bulk_dge(
  object,
  group_col = "case_control",
  .verbose = FALSE
)
object <- calculate_pca_bulk_dge(object, no_hvg_genes = 500L)
object <- batch_correction_bulk_dge(
  object,
  contrast_column = "case_control",
  batch_col = "batch",
  no_hvg_genes = 500L
)
get_outputs(object)$normalised_counts_corrected[1:3, 1:3]
#>        sample_1 sample_10 sample_100
#> gene_1 9.486898  8.890003   8.914310
#> gene_2 3.838447  2.857655   5.938676
#> gene_3 6.054164  5.771220   9.178742
```
