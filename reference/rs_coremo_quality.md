# Helper function to assess CoReMo cluster quality

This function assesses the quality of the clusters with a given cut `k`.
Returns the median R2 (cor^2) and the median absolute deviation (MAD) of
the clusters. Large clusters (≥1000) are subsampled to a random set of
1000 genes.

## Usage

``` r
rs_coremo_quality(cluster_genes, cor_mat, row_names, seed)
```

## Arguments

- cluster_genes:

  A list. Contains the cluster and their respective genes.

- cor_mat:

  Numerical matrix. Contains the correlation coefficients.

- row_names:

  String vector. The row names (or column names) of the correlation
  matrix.

- seed:

  Integer. Random seed for the sub sampling of genes.

## Value

A list containing:

- r2med - median R2 of the cluster.

- r2mad - median absolute deviation of the R2 in the cluster.

- size - size of the cluster.
