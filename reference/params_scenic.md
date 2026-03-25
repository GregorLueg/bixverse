# Constructor for SCENIC parameters

Constructor for SCENIC parameters

## Usage

``` r
params_scenic(
  min_counts = 50L,
  min_cells = 0.03,
  learner_type = "randomforest",
  gene_batch_strategy = "correlated",
  gene_batch_size = NULL,
  n_pcs = 50L,
  n_subsample = 100000L,
  learner_params = list()
)
```

## Arguments

- min_counts:

  Integer. Minimum total counts a gene needs to be included in the
  analysis. Defaults to `50L`.

- min_cells:

  Numeric. Minimum proportion of cells (between 0 and 1) that must
  express a gene for it to be considered. Defaults to `0.03`.

- learner_type:

  Character. Regression learner to use. One of `"randomforest"`,
  `"extratrees"`, or `"grnboost2"`. Defaults to `"randomforest"`.

- gene_batch_strategy:

  Character. Strategy for grouping target genes into batches. One of
  `"random"` or `"correlated"`. Only used for `"randomforest"` and
  `"extratrees"` learners; ignored for `"grnboost2"`. Defaults to
  `"correlated"`.

- gene_batch_size:

  Optional integer. Number of genes per batch. If `NULL` (default), the
  batch size is determined automatically. Ignored for `"grnboost2"`.

- n_pcs:

  Integer. Number of PCs to use for the correlated gene batch strategy.
  Defaults to `50L`.

- n_subsample:

  Integer. Cell subsampling threshold for the correlated gene batch
  strategy. If the number of cells meets or exceeds this value,
  `n_subsample` cells are randomly selected prior to running randomised
  SVD. Defaults to `100000L`.

- learner_params:

  List. Optional overrides for the regression learner parameters. For
  `"randomforest"`, see
  [`params_scenic_random_forest_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_scenic_random_forest_defaults.md).
  For `"extratrees"`, see
  [`params_scenic_extra_trees_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_scenic_extra_trees_defaults.md).
  For `"grnboost2"`, see
  [`params_scenic_gradient_boosting_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_scenic_gradient_boosting_defaults.md).

## Value

A named flat list with all SCENIC parameters.
