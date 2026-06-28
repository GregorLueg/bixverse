# Run scDblFinder doublet detection on a SingleCells object

Cluster-aware doublet detection using engineered features and a
gradient-boosted classifier. See Germain et al., F1000Research, 2022.

## Usage

``` r
scdblfinder_sc(
  object,
  scdblfinder_params = params_scdblfinder(),
  return_features = FALSE,
  cells_to_use = NULL,
  group_by = NULL,
  streaming = NULL,
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` class.

- scdblfinder_params:

  List. Parameters from
  [`params_scdblfinder()`](https://gregorlueg.github.io/bixverse/reference/params_scdblfinder.md).

- return_features:

  Boolean. Shall the features used to train the classifier be returned.

- cells_to_use:

  Optional string. Names of the cells to use for the run of the boosted
  doublet detection. Useful when you wish to run doublet detection on
  individual batches within your data. The object returned will be
  specifically using these cells.

- group_by:

  Optional grouping variable. Useful if you want to run the method on a
  per-sample basis.

- streaming:

  Optional Boolean. Shall the data be streamed in. Useful for larger
  data sets where you wish to avoid loading in the whole data. If
  `NULL`, will automatically detect.

- seed:

  Integer. Seed for reproducibility.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

An S3 object of class `ScDblFinderRes` containing:

- predicted_doublets:

  Logical vector of doublet calls.

- doublet_score:

  Numeric vector of classifier probabilities.

- cxds_scores:

  Numeric vector of the cxds scores.

- weighted:

  Numeric vector of the weighted scores.

- threshold:

  The threshold used for calling.

- cluster_labels:

  Integer vector of final cluster assignments.

- detected_doublet_rate:

  Fraction of cells called as doublets.

with `cell_indices` stored as an attribute.
