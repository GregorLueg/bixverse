# Row-bind the count matrices of meta cell objects

The inputs are already CSR, so a row-bind is concatenation of their
`@p`, `@j` and `@x` slots. No COO round-trip and no
`Matrix::sparseMatrix(i, j, x)` call, which would sort and copy the
whole index structure once per assay. The output vectors are allocated
once off a counting pass, so peak memory is the result plus one input's
scratch.

## Usage

``` r
.merge_mc_counts(inputs, gene_maps, n_genes)
```

## Arguments

- inputs:

  List of `MetaCells` objects.

- gene_maps:

  List of integer vectors mapping each input's gene positions onto the
  target gene space, 1-indexed. `NA` marks genes to drop.

- n_genes:

  Integer. Number of genes in the target gene space.

## Value

A list in the shape
[`get_meta_cell_matrices()`](https://gregorlueg.github.io/bixverse/reference/get_meta_cell_matrices.md)
consumes:

- indptr - Integer vector. 0-indexed CSR row pointer.

- indices - Integer vector. 0-indexed column indices.

- raw_counts - Numeric vector. Raw counts.

- norm_counts - Numeric vector. Normalised counts.

- nrow - Integer. Number of merged meta cells.

- ncol - Integer. Number of genes.
