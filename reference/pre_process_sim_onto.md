# Pre-process data for subsequent ontology similarity

This function calculates needed information for semantic similiary
calculations

## Usage

``` r
pre_process_sim_onto(object, .verbose = TRUE)
```

## Arguments

- object:

  `OntologySim class`. See
  [`OntologySim()`](https://gregorlueg.github.io/bixverse/reference/OntologySim.md).

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

The class with added pre-processed data for semantic similarities to the
properties.

## Examples

``` r
# ancestors, descendants and information content added to the class
onto <- data.table::data.table(
  parent = c("a", "b", "b", "b", "c"),
  child = c("b", "c", "d", "e", "f"),
  type = c("part_of", "part_of", "part_of", "is_a", "is_a")
)
onto_obj <- OntologySim(onto, .verbose = FALSE)
onto_obj <- pre_process_sim_onto(onto_obj, .verbose = FALSE)
names(S7::prop(onto_obj, "outputs"))
#> [1] "ancestors"           "descendants"         "information_content"
```
