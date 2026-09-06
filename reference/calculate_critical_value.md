# Calculates the critical value

This function calculates the critical value for a given ontology
similarity matrix.

## Usage

``` r
calculate_critical_value(x, alpha, permutations = 100000L, seed = 10101L)
```

## Arguments

- x:

  Numerical matrix or `OntologySim` class, see
  [`OntologySim()`](https://gregorlueg.github.io/bixverse/reference/OntologySim.md).
  This function tends to be slower on matrices compared to the
  `OntologySim` class.

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

## Examples

``` r
# critical value of a Wang similarity matrix at alpha 0.1
onto <- data.table::data.table(
  parent = c("a", "b", "b", "b", "c"),
  child = c("b", "c", "d", "e", "f"),
  type = c("part_of", "part_of", "part_of", "is_a", "is_a")
)
weights <- c(part_of = 0.8, is_a = 0.6)
sim_mat <- calculate_wang_sim_mat(onto, weights = weights)
calculate_critical_value(sim_mat, alpha = 0.1)
#> [1] 0.7641509
```
