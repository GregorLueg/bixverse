# Set similarities over two list

**\[experimental\]** This function calculates the Jaccard similarity or
overlap coefficient between all elements of two lists.

## Usage

``` r
rs_set_similarity_list2(s_1_list, s_2_list, overlap_coefficient)
```

## Arguments

- s_1_list:

  R list. The first list of string elements you want to compare against.

- s_2_list:

  R list. The second list of string elements you want to compare
  against.

- overlap_coefficient:

  Boolean. Use the overlap coefficient instead of the Jaccard
  similarity.

## Value

A matrix of the similarities between the elements. The rows represent
`s_1_list` and the columns `s_2_list`.
