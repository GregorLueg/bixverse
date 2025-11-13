# Filter the term similarities for a specific critical value

This function takes the similarity values as the upper triangle, the
row/column names and filtering the values down based on the threshold.

## Usage

``` r
rs_filter_onto_sim(sim_vals, names, threshold)
```

## Arguments

- sim_vals:

  Numerical vector. The upper triangle of the similarity matrix as a
  flattened vector.

- names:

  String vector. The row/col names of the similarity matrix.

- threshold:

  Float. The filtering threshold.

## Value

A list with:

- t1 - name of term 1.

- t2 - name of term 2.

- sim - the similarity between the two terms.
