# Detect Doublets via BoostClassifier (in Rust)

Detect Doublets via BoostClassifier (in Rust)

## Usage

``` r
rs_sc_doublet_detection(
  f_path_gene,
  f_path_cell,
  cells_to_keep,
  boost_params,
  seed,
  streaming,
  verbose
)
```

## Arguments

- f_path_gene:

  String. Path to the `counts_genes.bin` file.

- f_path_cell:

  String. Path to the `counts_cells.bin` file.

- cells_to_keep:

  Integer vector. The indices (0-indexed!) of the cells to include in
  this analysis.

- boost_params:

  List. Parameter list, see [`params_boost()`](params_boost.md).

- seed:

  Integer. Seed for reproducibility purposes.

- streaming:

  Boolean. Shall the data be streamed for the HVG calculations.

- verbose:

  Boolean. Controls verbosity

## Value

A list with

- predicted_doublets - Boolean vector indicating which observed cells
  predicted as doublets (TRUE = doublet, FALSE = singlet).

- doublet_scores_obs - Numerical vector with the likelihood of being a
  doublet for the observed cells.

- voting_avg - Voting average across the different iterations.
