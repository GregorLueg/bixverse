# OntologySim class

This class is used to store any ontology and apply different methods to
it. Currently implemented are different types of term similarities,
i.e., based on semantic similarities and the Wang similarity.

## Usage

``` r
OntologySim(parent_child_dt, .verbose = TRUE)
```

## Arguments

- parent_child_dt:

  A data.table that contains the ontological information in terms of
  parent child relationships. Need to contain the `c("parent", "child")`
  columns.

- .verbose:

  Boolean. Controls the verbosity of the class

## Value

Returns the class for subsequent usage.

## Properties

- edge_dt:

  data.table. Contains the parent-child relationships. (For Wang
  similarity also the relationship type.)

- outputs:

  List. Contains various intermediary results used for some methods.

- sim_mat:

  List. Contains the potentially calculated similarity matrix in form of
  an R6 class. Getters to access the data are provided.

- params:

  A (nested) list that will store all the parameters of the applied
  function.

- final_results:

  Final results stored in the class.

## Examples

``` r
# wrap a toy parent-child ontology into the class
onto <- data.table::data.table(
  parent = c("a", "b", "b", "b", "c"),
  child = c("b", "c", "d", "e", "f"),
  type = c("part_of", "part_of", "part_of", "is_a", "is_a")
)
OntologySim(onto, .verbose = FALSE)
#> OntologySim class:
#>  Size of the ontology: 5.
#>  Semantic similarities calculated: No.
```
