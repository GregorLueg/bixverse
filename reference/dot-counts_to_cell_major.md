# Reinterpret a gene-major count matrix as a cell-major one

Seurat and `SingleCellExperiment` both store counts as a `dgCMatrix` of
genes x cells. Everything on this side wants cells x genes. Those two
are the same bytes: a CSC matrix over genes x cells holds one pointer
per cell and gene indices, and so does a CSR matrix over cells x genes.
So this relabels the slots rather than transposing, and costs nothing
regardless of how many non-zeros there are.

## Usage

``` r
.counts_to_cell_major(counts)
```

## Arguments

- counts:

  `dgCMatrix`. Counts of genes x cells.

## Value

The same data as a `dgRMatrix` of cells x genes.
