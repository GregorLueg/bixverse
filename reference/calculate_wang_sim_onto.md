# Calculate the Wang similarity for an ontology.

This function calculates the Wang similarity for the whole ontology and
adds it to the class in a memory-efficient format for subsequent usage.

## Usage

``` r
calculate_wang_sim_onto(object, weights, .verbose = TRUE)
```

## Arguments

- object:

  `OntologySim` class. See
  [`OntologySim()`](https://gregorlueg.github.io/bixverse/reference/OntologySim.md).

- weights:

  Named numeric. The relationship of type to weight for this specific
  edge. For example `c("part_of" = 0.8, "is_a" = 0.6)`.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

The class with added semantic similarities to the properties.

## Examples

``` r
# Wang similarities for the whole ontology stored in the class
onto <- data.table::data.table(
  parent = c("a", "b", "b", "b", "c"),
  child = c("b", "c", "d", "e", "f"),
  type = c("part_of", "part_of", "part_of", "is_a", "is_a")
)
onto_obj <- OntologySim(onto, .verbose = FALSE)
calculate_wang_sim_onto(
  onto_obj,
  weights = c(part_of = 0.8, is_a = 0.6),
  .verbose = FALSE
)
#> OntologySim class:
#>  Size of the ontology: 5.
#>  Semantic similarities calculated: Yes.
```
