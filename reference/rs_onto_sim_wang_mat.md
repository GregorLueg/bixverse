# Calculate the Wang similarity matrix for an ontology

This function calculates the Wang similarity matrix for a given
ontology.

## Usage

``` r
rs_onto_sim_wang_mat(parents, children, w, flat_matrix)
```

## Arguments

- parents:

  String vector. The names of the parents.

- children:

  String vector. The names of the childs. The length of `parents` needs
  to be equal to `children`.

- w:

  Numerics. The weights between the parents and children. Need to be
  values between 0 and 1.

- flat_matrix:

  Boolean. Shall only the upper triangle be returned.

## Value

A list with:

- sim_mat - the Wang similarity matrix (flat or as matrix.)

- names - the row and column names for the calculated matrix.
