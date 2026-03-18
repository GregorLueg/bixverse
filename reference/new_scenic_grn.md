# Constructor for SCENIC GRN results

Stores the TF-gene importance matrix and associated metadata from a
SCENIC GRN inference run. Intended as input to downstream regulon
generation functions.

## Usage

``` r
new_scenic_grn(importance_matrix, gene_ids, tf_ids, params)
```

## Arguments

- importance_matrix:

  Matrix. Importance matrix of shape `(n_genes, n_tfs)` with gene
  identifiers as rownames and TF identifiers as colnames.

- gene_ids:

  Character vector. Gene identifiers corresponding to rows.

- tf_ids:

  Character vector. TF identifiers corresponding to columns.

- params:

  List. The full SCENIC parameters used for the run.

## Value

An object of class `ScenicGrn`.
