# Calculate the Resnik or Lin semantic similarity matrix

This function calculates the semantic similarities based on Resnik or
Lin similarity for a given ontology. Has also the option to calculate a
combined version of the two. Returns the full matrix of the data.

## Usage

``` r
calculate_semantic_sim_mat(similarity_type, ancestor_list, ic_list)
```

## Arguments

- similarity_type:

  String. One of `c("resnik", "lin", "combined")`. The Resnik similarity
  is automatically rescaled between 0 and 1 (dividing the values by
  maximum information content observed).

- ancestor_list:

  List. Names being the term and the elements in the list the names of
  the ancestors, see
  [`get_ontology_ancestry()`](get_ontology_ancestry.md).

- ic_list:

  List. The names being the term and the elements the information
  content of this given term. Needs to be a single float! See
  [`calculate_information_content()`](calculate_information_content.md).

## Value

The symmetric similarity matrix for the specified data from the
ontology.
