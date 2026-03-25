# Run SCENIC GRN inference

Runs SCENIC GRN inference on the provided genes using the specified
transcription factors as predictors. Returns a `ScenicGrn` object
containing the TF-gene importance matrix for further processing.

## Usage

``` r
scenic_grn_sc(
  object,
  tf_ids,
  scenic_params = params_scenic(),
  genes_to_take = NULL,
  cells_to_take = NULL,
  streaming = FALSE,
  random_seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` class.

- tf_ids:

  Character vector. Transcription factor gene identifiers to use as
  predictors. Must be a subset of gene identifiers present in the
  object.

- scenic_params:

  List. SCENIC parameters, see
  [`params_scenic()`](https://gregorlueg.github.io/bixverse/reference/params_scenic.md).

- genes_to_take:

  Optional character vector. Target gene identifiers. If `NULL`, genes
  are selected automatically via
  [`scenic_gene_filter_sc()`](https://gregorlueg.github.io/bixverse/reference/scenic_gene_filter_sc.md)
  using the `min_counts` and `min_cells` thresholds in `scenic_params`.

- cells_to_take:

  Optional string vector. Cell identifiers to restrict to. If `NULL`,
  defaults to all filtered cells in the class.

- streaming:

  Boolean. Whether to use the streaming implementation to bound memory
  usage. Useful for large datasets. Defaults to `FALSE`.

- random_seed:

  Integer. Used for reproducibility. Defaults to `42L`.

- .verbose:

  Boolean. Controls verbosity. Defaults to `TRUE`.

## Value

A `ScenicGrn` object.

## Details

TF identifiers that are not found in the object's gene list are silently
dropped with a warning indicating how many were removed. TF indices are
intersected with the target gene indices so that TFs not passing the
gene filter are excluded from the predictor set but remain as potential
targets if present in `genes_to_take`. You have the option to generate
the TF-gene importance values with three distinct methods. For the
`random_forest` and the `extratrees` version, a batching strategy is
applied in the default settings. Correlated genes are identified and
clustered together via k-means clustering on the feature loadings of the
PCA. These are then divided into batches of `gene_batch_size` and the
regression learners are leveraging multi-target regression to fit all
genes in the batch in one go. This massively accelerates the algorithm
and the importance values per gene-TF pair are calculated then
individually. Due to the batching by similar gene, the signal dilution
is limited. If you wish to run the traditional approach, you can set
gene_batch_size to `1L` or use the `grnboost2` learner that can only fit
one gene at a given time.
