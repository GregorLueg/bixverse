# Get TF gene correlations for MetaCells

Get TF gene correlations for MetaCells

## Usage

``` r
.tf_gene_cor_mc(tf_to_gene, object, spearman)
```

## Arguments

- tf_to_gene:

  data.table. The TF to gene data.table

- object:

  `MetaCells` class.

- spearman:

  Boolean. Shall the Spearman correlation be used.

## Value

Adds a `pairwise_cor` column to the tf_to_gene data.table
