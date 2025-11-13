# Set similarities over one list

This function calculates the set similarity via Jaccard or overlap
coefficient across all permutations of one list.

## Usage

``` r
rs_set_similarity_list(list, overlap_coefficient)
```

## Arguments

- list:

  A named R list.

- overlap_coefficient:

  Boolean. Use the overlap coefficient instead of the Jaccard similarity
  be calculated.

## Value

A list with the following items:

- from - Name of element i

- to - Name of element j

- sim - Similarity between the two elements
