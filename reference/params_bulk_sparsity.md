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
  strategy = c("seq_depth"),
  target_library_size = 20000,
  capture_efficiency_sigma = 0.5,
  seed = 123L
)
```

## Arguments

- strategy:

  String. Which dropout strategy to apply. Currently only `"seq_depth"`.

- target_library_size:

  Float. Reference library size per sample.

- capture_efficiency_sigma:

  Float. Standard deviation of the LogNormal size-factor distribution.
  Larger values spread the library sizes further apart.

- seed:

  Integer. Seed for reproducibility purposes.

## Value

A list with the parameters for usage in subsequent functions.

## References

Zappia, et al., Genome Biol, 2017
