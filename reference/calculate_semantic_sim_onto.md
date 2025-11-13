# Calculate the Resnik or Lin semantic similarity for an ontology.

This function calculates the specified semantic similarities for the
whole ontology and adds it to the class.

## Usage

``` r
calculate_semantic_sim_onto(
  object,
  sim_type = c("resnik", "lin", "combined"),
  .verbose = TRUE
)
```

## Arguments

- object:

  `ontology class`. See [`ontology()`](ontology.md).

- sim_type:

  String. One of `c("resnik", "lin", "combined")`.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

The class with added semantic similarities to the properties.
