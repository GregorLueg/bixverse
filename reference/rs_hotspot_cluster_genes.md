# Cluster the genes by Z-score together

Cluster the genes by Z-score together

## Usage

``` r
rs_hotspot_cluster_genes(z_matrix, fdr_threshold, min_size)
```

## Arguments

- z_matrix:

  Numerical matrix representing the Z-scores.

- fdr_threshold:

  Float. The FDR thresholds in terms of the Z-scores.

- min_size:

  Integer. Minimum cluster size.

## Value

An assignment vector. NA indicates that the gene did not pass the
thresholds and has not been assigned.
