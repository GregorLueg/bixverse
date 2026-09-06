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

## Examples

``` r
# keep only term pairs above the permutation-derived critical value
onto <- data.table::data.table(
  parent = c("a", "b", "b", "b", "c"),
  child = c("b", "c", "d", "e", "f"),
  type = c("part_of", "part_of", "part_of", "is_a", "is_a")
)
onto_obj <- OntologySim(onto, .verbose = FALSE)
onto_obj <- calculate_wang_sim_onto(
  onto_obj,
  weights = c(part_of = 0.8, is_a = 0.6),
  .verbose = FALSE
)
onto_obj <- filter_similarities(onto_obj, alpha = 0.1, .verbose = FALSE)
head(get_results(onto_obj))
#>        t1     t2       sim
#>    <char> <char>     <num>
#> 1:      f      c 0.7960848
#> 2:      d      b 0.7641509
#> 3:      b      c 0.7641509
```
