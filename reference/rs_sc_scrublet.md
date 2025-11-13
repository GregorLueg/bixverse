# Scrublet Rust interface

Scrublet Rust interface

## Usage

``` r
rs_sc_scrublet(
  f_path_gene,
  f_path_cell,
  cells_to_keep,
  scrublet_params,
  seed,
  verbose,
  streaming,
  return_combined_pca,
  return_pairs
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

- scrublet_params:

  List. Parameter list, see [`params_scrublet()`](params_scrublet.md).

- seed:

  Integer. Seed for reproducibility purposes.

- verbose:

  Boolean. Controls verbosity

- streaming:

  Boolean. Shall the data be streamed for the HVG calculations.

- return_combined_pca:

  Boolean. Shall the generated PCA be returned.

- return_pairs:

  Boolean. Shall the parents of the simulated cells be returned.

## Value

A list with

- predicted_doublets - Boolean vector indicating which observed cells
  predicted as doublets (TRUE = doublet, FALSE = singlet).

- doublet_scores_obs - Numerical vector with the likelihood of being a
  doublet for the observed cells.

- doublet_scores_sim - Numerical vector with the likelihood of being a
  doublet for the simulated cells.

- doublet_errors_obs - Numerical vector with the standard errors of the
  scores for the observed cells.

- z_scores - Z-scores for the observed cells. Represents:
  `score - threshold / error`.

- threshold - Used threshold.

- detected_doublet_rate - Fraction of cells that are called as doublet.

- detectable_doublet_fraction - Fraction of simulated doublets with
  scores above the threshold.

- overall_doublet_rate - Estimated overall doublet rate. Should roughly
  match the expected doublet rate.

- pca - Optional PCA embeddings across the original cells and simulated
  doublets.

- pair_1 - Optional integer vector representing the first parent of the
  simulated doublets.

- pair_2 - Optional integer vector representing the second parent of the
  simulated doublets.
