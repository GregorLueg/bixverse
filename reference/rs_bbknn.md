# BBKNN implementation in Rust

This function implements the BBKNN algorithm from Polański, et al.

## Usage

``` r
rs_bbknn(embd, batch_labels, bbknn_params, seed, verbose)
```

## Arguments

- embd:

  Numerical matrix. The embedding matrix to use to generate the BBKNN
  parameters. Usually PCA. Rows represent cells.

- batch_labels:

  Integer vector. These represent to which batch a given cell belongs.

- bbknn_params:

  List. Contains all of the BBKNN parameters.

- seed:

  Integer. Seed for reproducibility purposes.

- verbose:

  Boolean. Controls verbosity of the function.

## Value

A list of two lists representing the sparse matrix representation of the
distances and the connectivities.

## References

Polański, et al., Bioinformatics, 2020
