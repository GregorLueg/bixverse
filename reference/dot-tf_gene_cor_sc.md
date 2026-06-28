# Get TF gene correlations for SingleCells

Get TF gene correlations for SingleCells

## Usage

``` r
.tf_gene_cor_sc(tf_to_gene, object, spearman)
```

## Arguments

- tf_to_gene:

  data.table. The TF to gene data.table

- object:

  `SingleCells` class.

- spearman:

  Boolean. Shall the Spearman correlation be used.

## Value

Adds a `pairwise_cor` column to the tf_to_gene data.table
