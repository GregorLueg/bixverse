# Wrapper function to generate bulk sparsification parameters

Parameters for
[`simulate_dropouts()`](https://gregorlueg.github.io/bixverse/reference/simulate_dropouts.md).
Dropout falls out of the library size rather than an explicit per-gene
dropout curve: a size factor
`s_j ~ LogNormal(0, capture_efficiency_sigma)` is drawn per sample,
giving a target library size of `target_library_size * s_j`, and each
gene is binomially thinned towards that target.

## Usage

``` r
params_bulk_sparsity(
  strategy = "seq_depth",
  target_library_size = 20000,
  capture_efficiency_sigma = 0.5,
  seed = 123L
)
```

## Arguments

- strategy:

  String. Which dropout strategy to apply. Currently only `"seq_depth"`.
  One of `"seq_depth"`. Defaults to `"seq_depth"`.

- target_library_size:

  Numeric. Reference library size per sample. Defaults to `20000.0`.

- capture_efficiency_sigma:

  Numeric. Standard deviation of the LogNormal size-factor distribution.
  Larger values spread the library sizes further apart. Defaults to
  `0.5`.

- seed:

  Integer. Seed for reproducibility purposes. Defaults to `123L`.

## Value

A named list with the following elements:

- strategy - String. Which dropout strategy to apply. Currently only
  `"seq_depth"`. One of `"seq_depth"`. Defaults to `"seq_depth"`.

- target_library_size - Numeric. Reference library size per sample.
  Defaults to `20000.0`.

- capture_efficiency_sigma - Numeric. Standard deviation of the
  LogNormal size-factor distribution. Larger values spread the library
  sizes further apart. Defaults to `0.5`.

- seed - Integer. Seed for reproducibility purposes. Defaults to `123L`.

## References

Zappia, et al., Genome Biol, 2017
