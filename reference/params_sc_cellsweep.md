# Default parameters for CellSweep denoising

Mirrors the CellSweep reference implementation's defaults. The
pseudocounts (`celltype_lambda`, `ambient_lambda`, `bulk_lambda`) are
given on the scale you see here and divided by the gene count
internally. Two of these are worth knowing about before you touch
anything else. `freeze_ambient_profile = TRUE` keeps the ambient profile
at its empty-droplet estimate, which is the recommended path and the
only one where `alpha_cap`, the repulsion terms and cell-type
reassignment are live. And `freeze_empties` only accepts `TRUE`: the
reference gives empty droplets a cell-type component they have no label
for, which indexes past the end of the profile matrix and wraps onto the
last cell type.

## Usage

``` r
params_sc_cellsweep(
  freeze_empties = TRUE,
  freeze_ambient_profile = TRUE,
  init_alpha = 0.9,
  init_beta = 0.1,
  alpha_cap = 0.9,
  repulsion_strength = 1e-04,
  max_frac_gene_repulsion = 0.2,
  celltype_lambda = 50,
  ambient_lambda = 50,
  bulk_lambda = 10,
  eps = 1e-12,
  log_eps = as.numeric("1e-300"),
  max_iter = 2000L,
  del0_ll_tol = 0.001,
  min_ll_tol = 1e-06,
  tol_p = 1e-04,
  tol_f = 1e-04,
  norm_from_rounded = FALSE,
  seed = 42L
)
```

## Arguments

- freeze_empties:

  Boolean. Keep the contamination fraction of empty droplets pinned
  at 1. Only `TRUE` is supported, see the description. Defaults to
  `TRUE`.

- freeze_ambient_profile:

  Boolean. Keep the ambient profile at its empty-droplet estimate rather
  than re-estimating it as a mixture over the cell-type profiles.
  Defaults to `TRUE`.

- init_alpha:

  Numeric. Starting ambient fraction for every real barcode. With
  `freeze_ambient_profile = TRUE` the final result barely depends on it,
  so it sits at `alpha_cap`. Defaults to `0.9`.

- init_beta:

  Numeric. Starting bulk contamination fraction. Set below `init_alpha`
  on purpose: bulk and ambient are not fully separable, so this biases
  unassignable contamination towards ambient. Defaults to `0.1`.

- alpha_cap:

  Numeric. Ceiling on the per-cell ambient fraction before the
  log-likelihood converges. Barcodes wanting to exceed it are excluded
  from the cell-type profile update and allowed to switch cell type.
  Defaults to `0.9`.

- repulsion_strength:

  Numeric. Strength of the repulsion pushing cell-type profiles away
  from the ambient profile. Scales with cluster mass, so it is inert on
  small data and only bites at realistic cell counts. Defaults to
  `1e-04`.

- max_frac_gene_repulsion:

  Numeric. Ceiling on the fraction of any single profile entry that
  repulsion may remove. Defaults to `0.2`.

- celltype_lambda:

  Numeric. Pseudocount smoothing the cell-type profile update. Higher
  values give smoother profiles. Defaults to `50.0`.

- ambient_lambda:

  Numeric. Pseudocount smoothing the ambient profile estimate. Defaults
  to `50.0`.

- bulk_lambda:

  Numeric. Pseudocount smoothing the bulk profile estimate. Defaults to
  `10.0`.

- eps:

  Numeric. Floor on denominators. Defaults to `1e-12`.

- log_eps:

  Numeric. Floor on the argument of `log`. Defaults to `1e-300`.

- max_iter:

  Integer. Hard cap on EM iterations. Defaults to `2000L`.

- del0_ll_tol:

  Numeric. Log-likelihood change, as a fraction of the first EM step's
  change, below which stage one ends and parameter convergence starts
  being checked. Defaults to `0.001`.

- min_ll_tol:

  Numeric. Floor on the adaptive tolerance, relative to the current
  log-likelihood. Stops `del0_ll_tol` chasing floating point noise.
  Defaults to `1e-06`.

- tol_p:

  Numeric. Convergence threshold on the maximum row-wise L1 change in
  the cell-type profiles. Defaults to `1e-04`.

- tol_f:

  Numeric. Convergence threshold on the change in the total
  contamination fraction. Defaults to `1e-04`.

- norm_from_rounded:

  Boolean. Derive the normalised layer from the integerised counts
  rather than the denoised floats. Consistent across the two layers at
  the cost of the sub-integer signal, which is where CellSweep is most
  informative. Defaults to `FALSE`.

- seed:

  Integer. Seed for the stochastic rounding of the denoised counts.
  Defaults to `42L`.

## Value

A named list with the following elements:

- freeze_empties - Boolean. Keep the contamination fraction of empty
  droplets pinned at 1. Only `TRUE` is supported, see the description.
  Defaults to `TRUE`.

- freeze_ambient_profile - Boolean. Keep the ambient profile at its
  empty-droplet estimate rather than re-estimating it as a mixture over
  the cell-type profiles. Defaults to `TRUE`.

- init_alpha - Numeric. Starting ambient fraction for every real
  barcode. With `freeze_ambient_profile = TRUE` the final result barely
  depends on it, so it sits at `alpha_cap`. Defaults to `0.9`.

- init_beta - Numeric. Starting bulk contamination fraction. Set below
  `init_alpha` on purpose: bulk and ambient are not fully separable, so
  this biases unassignable contamination towards ambient. Defaults to
  `0.1`.

- alpha_cap - Numeric. Ceiling on the per-cell ambient fraction before
  the log-likelihood converges. Barcodes wanting to exceed it are
  excluded from the cell-type profile update and allowed to switch cell
  type. Defaults to `0.9`.

- repulsion_strength - Numeric. Strength of the repulsion pushing
  cell-type profiles away from the ambient profile. Scales with cluster
  mass, so it is inert on small data and only bites at realistic cell
  counts. Defaults to `1e-04`.

- max_frac_gene_repulsion - Numeric. Ceiling on the fraction of any
  single profile entry that repulsion may remove. Defaults to `0.2`.

- celltype_lambda - Numeric. Pseudocount smoothing the cell-type profile
  update. Higher values give smoother profiles. Defaults to `50.0`.

- ambient_lambda - Numeric. Pseudocount smoothing the ambient profile
  estimate. Defaults to `50.0`.

- bulk_lambda - Numeric. Pseudocount smoothing the bulk profile
  estimate. Defaults to `10.0`.

- eps - Numeric. Floor on denominators. Defaults to `1e-12`.

- log_eps - Numeric. Floor on the argument of `log`. Defaults to
  `1e-300`.

- max_iter - Integer. Hard cap on EM iterations. Defaults to `2000L`.

- del0_ll_tol - Numeric. Log-likelihood change, as a fraction of the
  first EM step's change, below which stage one ends and parameter
  convergence starts being checked. Defaults to `0.001`.

- min_ll_tol - Numeric. Floor on the adaptive tolerance, relative to the
  current log-likelihood. Stops `del0_ll_tol` chasing floating point
  noise. Defaults to `1e-06`.

- tol_p - Numeric. Convergence threshold on the maximum row-wise L1
  change in the cell-type profiles. Defaults to `1e-04`.

- tol_f - Numeric. Convergence threshold on the change in the total
  contamination fraction. Defaults to `1e-04`.

- norm_from_rounded - Boolean. Derive the normalised layer from the
  integerised counts rather than the denoised floats. Consistent across
  the two layers at the cost of the sub-integer signal, which is where
  CellSweep is most informative. Defaults to `FALSE`.

- seed - Integer. Seed for the stochastic rounding of the denoised
  counts. Defaults to `42L`.
