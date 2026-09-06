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

  `OntologySim` class. See
  [`OntologySim()`](https://gregorlueg.github.io/bixverse/reference/OntologySim.md).

- sim_type:

  String. One of `c("resnik", "lin", "combined")`.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

The class with added semantic similarities to the properties.

## Examples

``` r
# Resnik similarities for the whole ontology stored in the class
onto <- data.table::data.table(
  parent = c("a", "b", "b", "b", "c"),
  child = c("b", "c", "d", "e", "f"),
  type = c("part_of", "part_of", "part_of", "is_a", "is_a")
)
onto_obj <- pre_process_sim_onto(
  OntologySim(onto, .verbose = FALSE),
  .verbose = FALSE
)
calculate_semantic_sim_onto(onto_obj, sim_type = "resnik", .verbose = FALSE)
#> OntologySim class:
#>  Size of the ontology: 5.
#>  Semantic similarities calculated: Yes.
```
