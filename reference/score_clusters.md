# Score clusters based on ScType

Score clusters based on ScType

## Usage

``` r
score_clusters(x, cluster_labels)

# S3 method for class 'ScTypeResults'
score_clusters(x, cluster_labels)
```

## Arguments

- x:

  `ScTypeResults` object.

- cluster_labels:

  Integer vector. Cluster assignment, of length of the scored cells.

## Value

A `data.table` with cluster_id, cell_type, scores and n_cells.
