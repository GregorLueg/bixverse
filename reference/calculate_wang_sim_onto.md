# Calculate the Wang similarity for an ontology.

This function calculates the Wang similarity for the whole ontology and
adds it to the class in a memory-efficient format for subsequent usage.

## Usage

``` r
calculate_wang_sim_onto(object, weights, .verbose = TRUE)
```

## Arguments

- object:

  `ontology class`. See [`ontology()`](ontology.md).

- weights:

  Named numeric. The relationship of type to weight for this specific
  edge. For example `c("part_of" = 0.8, "is_a" = 0.6)`.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

The class with added semantic similarities to the properties.
