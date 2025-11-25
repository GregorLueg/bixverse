# Cell cycle genes

(Human) cell cycle genes for scoring in single cell, please see
[`module_scores_sc()`](module_scores_sc.md) and Tirosh et al, Science
(2016). These genes are based on the updated version in Seurat.

## Usage

``` r
cell_cycle_genes
```

## Format

### `cell_cycle_genes`

A data.table

- hgnc_symbol:

  The (HGNC) gene symbol

- ensembl_gene_id:

  The corresponding ensembl identifiers.

- set:

  To which cell cycle phase the gene belongs
