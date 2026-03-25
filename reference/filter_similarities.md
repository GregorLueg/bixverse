# Filter the calculated similarities

This function calculates the critical value, see
[`calculate_critical_value()`](https://gregorlueg.github.io/bixverse/reference/calculate_critical_value.md)
and filters subsequently all the term pairs to the ones with a value ≥
critical value.

## Usage

``` r
filter_similarities(
  object,
  alpha,
  permutations = 100000L,
  seed = 10101L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `OntologySim` class. See
  [`OntologySim()`](https://gregorlueg.github.io/bixverse/reference/OntologySim.md).

- alpha:

  Float. The alpha value. For example, 0.001 would mean that the
  critical value is smaller than 0.1 percentile of the random
  permutations.

- permutations:

  Number of random permutations.

- seed:

  Integer. For reproducibility purposes.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

The class with filtered results added to the respective slot.
