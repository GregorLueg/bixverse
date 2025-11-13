# Doublet detection with Scrublet

This function implements the doublet detection from Scrublet, see
Wolock, et al. Briefly, arteficial doublets are being generated from the
data via random combination of initial cells. Highly variable genes
(HVG) are being identified and the observed cells are being projected on
a PCA space. Subsequently, the simulated doublets are being projected on
the same PCA space given the same HVGs; kNN graphs are being generated
and a kNN classifier is used to assign a probability that a given cell
in the original data is a doublet. For more details, please check the
publication.

## Usage

``` r
scrublet_sc(
  object,
  scrublet_params = params_scrublet(),
  seed = 42L,
  streaming = FALSE,
  return_combined_pca = FALSE,
  return_pairs = FALSE,
  .verbose = TRUE
)
```

## Arguments

- object:

  `single_cell_exp` class.

- scrublet_params:

  A list with the final scrublet parameters, see
  [`params_scrublet()`](params_scrublet.md) for full details.

- seed:

  Integer. Random seed.

- streaming:

  Boolean. Shall streaming be used during the HVG calculations. Slower,
  but less memory usage.

- return_combined_pca:

  Boolean. Shall the PCA of the observed cells and simulated doublets be
  returned.

- return_pairs:

  Boolean. Shall the pairs be returned.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

A `scrublet_res` class that has with the following items:

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

- pair_1 - Optional index of the parent cell 1 of the simulated
  doublets.

- pair_2 - Optional index of the parent cell 2 of the simulated
  doublets.

## References

Wollock, et al., Cell Syst, 2020
