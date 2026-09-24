# Wrapper function for analytic Pearson residual parameters

The closed-form alternative to scTransform: one shared dispersion
instead of a fitted model per gene. Much cheaper, and on most data sets
it ranks genes about as well.

## Usage

``` r
params_sc_apr(theta = 100, min_cells = 5L, clip_min = NULL, clip_max = NULL)
```

## Arguments

- theta:

  Numeric. The shared negative binomial dispersion. `Inf` gives the
  Poisson limit. Defaults to `100.0`.

- min_cells:

  Integer. Minimum number of cells a gene must be detected in to be
  retained. `0L` keeps everything. Defaults to `5L`.

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

- theta - Numeric. The shared negative binomial dispersion. `Inf` gives
  the Poisson limit. Defaults to `100.0`.

- min_cells - Integer. Minimum number of cells a gene must be detected
  in to be retained. `0L` keeps everything. Defaults to `5L`.

- clip_min - Numeric or `NULL`. Lower residual clipping bound. `NULL`
  uses `-sqrt(n_cells)`. Must be given together with `clip_max`.
  Defaults to `NULL`.

- clip_max - Numeric or `NULL`. Upper residual clipping bound. `NULL`
  uses `sqrt(n_cells)`. Must be given together with `clip_min`. Defaults
  to `NULL`.

## References

Lause, Berens and Kobak, Genome Biology, 2021.
