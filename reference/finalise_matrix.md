# Finalise the count matrix

Finalise the count matrix

## Usage

``` r
finalise_matrix(matrix, return_format, cell_indices, gene_indices, sc_map)
```

## Arguments

- matrix:

  The sparse matrix to finalise

- return_format:

  String. One of `c("cell", "gene")`.

- cell_indices:

  Optional integer. The cell indices to return.

- gene_indices:

  Optional integer. The gene indices to return.

- sc_map:

  A `sc_mapper` class. Contains various mapping information.

## Value

The finalised matrix.
