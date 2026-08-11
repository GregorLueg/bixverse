# Wrapper function for MAGIC imputation parameters

Parameters controlling the MAGIC imputation run by
[`run_magic_sc()`](https://gregorlueg.github.io/bixverse/reference/run_magic_sc.md).
The defaults mirror the reference implementation.

## Usage

``` r
params_sc_magic(
  n_steps = 3L,
  clip_threshold = 0.01,
  gene_batch_size = 1000L,
  layer = c("norm", "raw"),
  allow_large = FALSE
)
```

## Arguments

- n_steps:

  Integer. Diffusion steps applied to the counts. Defaults to `3L`. Zero
  is legal and hands back the un-imputed values, which is a cheap way to
  compare the two.

- clip_threshold:

  Numeric. Imputed values below this are zeroed after the last step.
  Defaults to `0.01`.

- gene_batch_size:

  Integer. Genes streamed off the binary store per block. Bounds the
  scratch memory and is clamped to the number of requested genes.
  Defaults to `1000L`.

- layer:

  String. One of `c("norm", "raw")`. Which stored layer to impute. The
  operator preserves per-cell mass, so imputed values sit on the scale
  of whatever went in: imputing raw counts and imputing log-normalised
  counts are different operations rather than the same one rescaled.
  Defaults to `"norm"`.

- allow_large:

  Boolean. Skip the output size guard. The dense output is capped at 1e9
  elements, i.e. 4 GB of `f32`. Defaults to `FALSE`.

## Value

A named flat list with all MAGIC parameters.

## References

van Dijk, et al., Cell, 2018.
