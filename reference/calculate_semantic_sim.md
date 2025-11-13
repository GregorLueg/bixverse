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
  [`get_ontology_ancestry()`](get_ontology_ancestry.md).

- ic_list:

  List. The names being the term and the elements the information
  content of this given term. Needs to be a single float! See
  [`calculate_information_content()`](calculate_information_content.md).

- add_self:

  Boolean. Shall self-similarities be added. Defaults to `FALSE`.

## Value

A data.table with the calculated similarities.
