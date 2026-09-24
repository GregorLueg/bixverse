# Wrapper function for scTransform (v2) parameters

Defaults are sctransform's own with `vst.flavor = "v2"` applied. Only
the step-1 fit scales with `n_genes` and `n_cells`: those bound the
subsample the negative binomial models are fitted on, and every later
pass streams gene by gene. Raising them costs fitting time, not memory.

## Usage

``` r
params_sc_sctransform(
  n_genes = 2000L,
  n_cells = 2000L,
  min_cells = 5L,
  bw_adjust = 3,
  gmean_eps = 1,
  outlier_th = 10,
  poisson_diff_theta = 0.001,
  clip_min = NULL,
  clip_max = NULL
)
```

## Arguments

- n_genes:

  Integer. Genes in the step-1 subsample. Defaults to `2000L`.

- n_cells:

  Integer. Cells in the step-1 subsample. Defaults to `2000L`.

- min_cells:

  Integer. Minimum number of cells a gene must be detected in to be
  modelled. Defaults to `5L`.

- bw_adjust:

  Numeric. Bandwidth multiplier for the kernel regression that
  regularises the parameters. Defaults to `3.0`.

- gmean_eps:

  Numeric. Offset in the geometric mean. Defaults to `1.0`.

- outlier_th:

  Numeric. Threshold, in median absolute deviations, past which a step-1
  fit is treated as an outlier. Defaults to `10.0`.

- poisson_diff_theta:

  Numeric. Below this, the fitted dispersion is taken as the Poisson
  limit. Defaults to `0.001`.

- clip_min:

  Numeric or `NULL`. Lower residual clipping bound. `NULL` uses
  `-sqrt(n_cells)`. Must be given together with `clip_max`. Defaults to
  `NULL`.

- clip_max:

  Numeric or `NULL`. Upper residual clipping bound. `NULL` uses
  `sqrt(n_cells)`. Must be given together with `clip_min`. Defaults to
  `NULL`.

## Value

A named list with the following elements:

- n_genes - Integer. Genes in the step-1 subsample. Defaults to `2000L`.

- n_cells - Integer. Cells in the step-1 subsample. Defaults to `2000L`.

- min_cells - Integer. Minimum number of cells a gene must be detected
  in to be modelled. Defaults to `5L`.

- bw_adjust - Numeric. Bandwidth multiplier for the kernel regression
  that regularises the parameters. Defaults to `3.0`.

- gmean_eps - Numeric. Offset in the geometric mean. Defaults to `1.0`.

- outlier_th - Numeric. Threshold, in median absolute deviations, past
  which a step-1 fit is treated as an outlier. Defaults to `10.0`.

- poisson_diff_theta - Numeric. Below this, the fitted dispersion is
  taken as the Poisson limit. Defaults to `0.001`.

- clip_min - Numeric or `NULL`. Lower residual clipping bound. `NULL`
  uses `-sqrt(n_cells)`. Must be given together with `clip_max`.
  Defaults to `NULL`.

- clip_max - Numeric or `NULL`. Upper residual clipping bound. `NULL`
  uses `sqrt(n_cells)`. Must be given together with `clip_min`. Defaults
  to `NULL`.

## References

Choudhary and Satija, Genome Biology, 2022.
