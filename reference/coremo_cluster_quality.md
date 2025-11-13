# Coremo: measures the quality of the clusters

Utility functions to give back the summary stats (median and mean R²)
for each of the identified clusters. Clusters over a size of 1000 genes
will be randomly sampled.

## Usage

``` r
coremo_cluster_quality(modules, cor_mat, random_seed = 10101L)
```

## Arguments

- modules:

  Named vector. The names reflect the gene of the associated module.

- cor_mat:

  Numeric matrix. The original correlation matrix.

- random_seed:

  Integer. Random seed to ensure consistency if sampling is used.

## Value

A data.table with the quality measures of the cluster.
