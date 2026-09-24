# Parameters for ligand to target influence computation

Parameters for ligand to target influence computation

## Usage

``` r
params_ligand_target(
  lr_sig_hub = 0,
  gr_hub = 0,
  ltf_cutoff = 0.99,
  damping_factor = 0.5,
  tol = 1e-06,
  max_iter = 1000L,
  topology_correction = FALSE,
  secondary_targets = FALSE
)
```

## Arguments

- lr_sig_hub:

  Numeric. Numeric in `[0, 1]`. Hub correction strength for the
  ligand-receptor / signalling layer. 0 disables correction. Defaults to
  `0.0`.

- gr_hub:

  Numeric. Numeric in `[0, 1]`. Hub correction strength for the gene
  regulatory layer. 0 disables correction. Defaults to `0.0`.

- ltf_cutoff:

  Numeric. Numeric in `[0, 1]`. Quantile cutoff applied to the
  intermediate ligand-to-TF matrix. Defaults to `0.99`.

- damping_factor:

  Numeric. Numeric in `[0, 1]`. PageRank-style damping factor. Defaults
  to `0.5`.

- tol:

  Numeric. Numeric \> 0. Convergence tolerance for the propagation step.
  Defaults to `1e-06`.

- max_iter:

  Integer. Integer \>= 1. Maximum iterations for the propagation step.
  Defaults to `1000L`.

- topology_correction:

  Boolean. Apply topology correction. Defaults to `FALSE`.

- secondary_targets:

  Boolean. Run a second round through targets. Defaults to `FALSE`.

## Value

A named list with the following elements:

- lr_sig_hub - Numeric. Numeric in `[0, 1]`. Hub correction strength for
  the ligand-receptor / signalling layer. 0 disables correction.
  Defaults to `0.0`.

- gr_hub - Numeric. Numeric in `[0, 1]`. Hub correction strength for the
  gene regulatory layer. 0 disables correction. Defaults to `0.0`.

- ltf_cutoff - Numeric. Numeric in `[0, 1]`. Quantile cutoff applied to
  the intermediate ligand-to-TF matrix. Defaults to `0.99`.

- damping_factor - Numeric. Numeric in `[0, 1]`. PageRank-style damping
  factor. Defaults to `0.5`.

- tol - Numeric. Numeric \> 0. Convergence tolerance for the propagation
  step. Defaults to `1e-06`.

- max_iter - Integer. Integer \>= 1. Maximum iterations for the
  propagation step. Defaults to `1000L`.

- topology_correction - Boolean. Apply topology correction. Defaults to
  `FALSE`.

- secondary_targets - Boolean. Run a second round through targets.
  Defaults to `FALSE`.
