# Resolve the genes DIALOGUE builds signatures from

Resolve the genes DIALOGUE builds signatures from

## Usage

``` r
.resolve_dialogue_genes(object, gene_ids)
```

## Arguments

- object:

  `SingleCells`, `SingleCellsSubset` or `MetaCells` class.

- gene_ids:

  Optional character. Genes to use. If `NULL`, uses the HVGs.

## Value

Integer vector with the 0-indexed gene positions.
