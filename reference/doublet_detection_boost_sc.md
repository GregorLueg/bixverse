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
  seed = 42L,
  streaming = FALSE,
  .verbose = TRUE
)
```

## Arguments

- object:

  `single_cell_exp` class.

- boost_params:

  A list with the final scrublet parameters, see
  [`params_boost()`](params_boost.md) for full details.

- seed:

  Integer. Random seed.

- streaming:

  Boolean. Shall streaming be used during the HVG calculations. Slower,
  but less memory usage.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

A `boost_res` class that has with the following items:

- predicted_doublets - Boolean vector indicating which observed cells
  predicted as doublets (TRUE = doublet, FALSE = singlet).

- doublet_scores_obs - Numerical vector with the likelihood of being a
  doublet for the observed cells.

- voting_avg - Numerical vector with the average voting score.
