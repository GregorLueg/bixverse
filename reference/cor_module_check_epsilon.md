# Iterate through different epsilon parameters

This functions iterates through a set of provided epsilons and checks
for each one to what extend the underlying affinity matrix will follow a
power law distribution.

## Usage

``` r
cor_module_check_epsilon(
  object,
  rbf_func = c("bump", "gaussian", "inverse_quadratic"),
  epsilons = c(0.25, seq(from = 0.5, to = 10, by = 0.5)),
  .verbose = TRUE
)
```

## Arguments

- object:

  The class, see [`bulk_coexp()`](bulk_coexp.md). You need to run
  [`cor_module_processing()`](cor_module_processing.md) before running
  this function.

- rbf_func:

  The type of RBF function you want to apply. A choice of
  `c('bump', 'gaussian', 'inverse_quadratic')`.

- epsilons:

  Vector of floats. The different epsilon parameters you would like to
  run.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

The class with added data to the properties for subsequent usage.
