# Calculate the Resnik or Lin semantic similarity

This function calculates the semantic similarities based on Resnik or
Lin similarity for a given ontology. Has also the option to calculate a
combined version of the two.

## Usage

``` r
calculate_semantic_sim(
  terms,
  similarity_type,
  ancestor_list,
  ic_list,
  add_self = FALSE
)
```

## Arguments

- terms:

  Vector of strings. The terms in the ontology you wish to screen.

- similarity_type:

  String. One of `c("resnik", "lin", "combined")`.

- ancestor_list:

  List. Names being the term and the elements in the list the names of
  the ancestors, see
  [`get_ontology_ancestry()`](https://gregorlueg.github.io/bixverse/reference/get_ontology_ancestry.md).

- ic_list:

  List. The names being the term and the elements the information
  content of this given term. Needs to be a single float! See
  [`calculate_information_content()`](https://gregorlueg.github.io/bixverse/reference/calculate_information_content.md).

- add_self:

  Boolean. Shall self-similarities be added. Defaults to `FALSE`.

## Value

A data.table with the calculated similarities.

## Examples

``` r
# Lin similarity for a subset of terms only
onto <- data.table::data.table(
  parent = c("a", "b", "b", "b", "c"),
  child = c("b", "c", "d", "e", "f")
)
ancestry <- get_ontology_ancestry(onto)
ic <- calculate_information_content(ancestry$descendants)
calculate_semantic_sim(
  terms = c("c", "d", "f"),
  similarity_type = "lin",
  ancestor_list = ancestry$ancestors,
  ic_list = ic
)
#>     term1  term2      sims
#>    <char> <char>     <num>
#> 1:      c      d 0.1261579
#> 2:      c      f 0.7601875
#> 3:      d      f 0.1017556
```
