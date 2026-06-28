# Score the individual clusters based on ScType

**\[experimental\]** This Rust function implements the cell type scoring
approach from Ianevski et al. (2022).

## Usage

``` r
rs_sc_type_cluster_assignment(sc_type_res, cluster_labels)
```

## Arguments

- sc_type_res:

  List. The ScType results.

- cluster_labels:

  Integer. Cluster assignment. Needs to be of length of scored cells.

## Value

A list with

- cluster_id - The cluster id/integer

- cell_type - String; the predicted cell type

- score - The final score for the clsuter.

- n_cells - The number of cells in the cluster.
