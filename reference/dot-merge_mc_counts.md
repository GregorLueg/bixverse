# Row-bind the count matrices of meta cell objects

Row-bind the count matrices of meta cell objects

## Usage

``` r
.merge_mc_counts(inputs, gene_maps, meta_cell_ids, gene_ids)
```

## Arguments

- inputs:

  List of `MetaCells` objects.

- gene_maps:

  List of integer vectors mapping each input's gene positions onto the
  target gene space. `NA` marks genes to drop.

- meta_cell_ids:

  Character vector with the merged meta cell identifiers.

- gene_ids:

  Character vector with the target gene identifiers.

## Value

A list with the merged `raw` and `norm` CSR matrices.
