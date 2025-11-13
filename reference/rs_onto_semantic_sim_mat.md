# Calculate the semantic similarity in an ontology

This function calculates the specified semantic similarity and returns
the full vector (only calculating the upper triangle) for the given
similarity.

## Usage

``` r
rs_onto_semantic_sim_mat(sim_type, ancestor_list, ic_list, flat_matrix)
```

## Arguments

- sim_type:

  String. Must be one of `c("resnik", "lin", "combined")`.

- ancestor_list:

  R list with names being the term and the elements in the list the
  names of the ancestors.

- ic_list:

  R list with the names being the term and the elements the information
  content of this given term. Needs to be a single float!

- flat_matrix:

  Boolean. Shall only the upper triangle be returned.

## Value

A list with:

- sim_mat - the semantic similarity matrix (flat or as matrix.)

- names - the row and column names for the calculated matrix.
