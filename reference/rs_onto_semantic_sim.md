# Calculate the semantic similarity in an ontology

This function calculates the specified semantic similarity and returns
the full vector (only calculating the upper triangle) for the given
similarity.

## Usage

``` r
rs_onto_semantic_sim(terms, sim_type, ancestor_list, ic_list)
```

## Arguments

- terms:

  Vector of strings. The terms in the ontology you wish to screen.

- sim_type:

  String. Must be one of `c("resnik", "lin", "combined")`.

- ancestor_list:

  R list with names being the term and the elements in the list the
  names of the ancestors.

- ic_list:

  R list with the names being the term and the elements the information
  content of this given term. Needs to be a single float!

## Value

A list with:

- term1 - name of the first term.

- term2 - name of the second term.

- sims - similarity between the two terms.
