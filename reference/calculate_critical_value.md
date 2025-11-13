# Calculates the critical value

This function calculates the critical value for a given ontology
similarity matrix.

## Usage

``` r
calculate_critical_value(x, alpha, permutations = 100000L, seed = 10101L)
```

## Arguments

- x:

  Numerical matrix or `ontology class`, see [`ontology()`](ontology.md).
  This function tends to be slower on matrices compared to
  `ontology class`.

- alpha:

  Float. The alpha value. For example, 0.001 would mean that the
  critical value is smaller than 0.1 percentile of the random
  permutations.

- permutations:

  Number of random permutations.

- seed:

  Integer. For reproducibility purposes

## Value

The critical value.
