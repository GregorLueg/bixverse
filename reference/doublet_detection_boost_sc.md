# Doublet detection with boosted doublet classification

This function implements the boosted doublet detection. It generates
through several iterations simulated doublets, generate kNN graphs, runs
Louvain clustering and assesses how often an observed cells clsuters
together with the simulated doublets.

## Usage

``` r
doublet_detection_boost_sc(
  object,
  boost_params = params_boost(),
  cells_to_use = NULL,
  group_by = NULL,
  seed = 42L,
  streaming = NULL,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` class.

- boost_params:

  A list with the final scrublet parameters, see
  [`params_boost()`](https://gregorlueg.github.io/bixverse/reference/params_boost.md)
  for full details.

- cells_to_use:

  Optional string. Names of the cells to use for the run of the boosted
  doublet detection. Useful when you wish to run doublet detection on
  individual batches within your data. The object returned will be
  specifically using these cells.

- group_by:

  Optional grouping variable. Useful if you want to run the method on a
  per-sample basis.

- seed:

  Integer. Random seed.

- streaming:

  Optional Boolean. Shall the data be streamed in. Useful for larger
  data sets where you wish to avoid loading in the whole data. If
  `NULL`, will automatically detect.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

A `boost_res` class that has with the following items:

- predicted_doublets - Boolean vector indicating which observed cells
  predicted as doublets (TRUE = doublet, FALSE = singlet).

- doublet_scores_obs - Numerical vector with the likelihood of being a
  doublet for the observed cells.

- voting_avg - Numerical vector with the average voting score.

## Examples

``` r
# boosted doublet detection over five iterations
sc <- demo_single_cells(prepped = FALSE)
doublet_detection_boost_sc(
  sc,
  boost_params = params_boost(
    hvg = list(min_gene_var_pctl = 0.0),
    pca = list(no_pcs = 10L),
    n_iters = 5L
  ),
  .verbose = FALSE
)
#> BoostRes: 500 cells, 1 doublets (0.2%)
#>   Score range: [0.0125, 0.8289]

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
