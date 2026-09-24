# Score the individual clusters based on ScType

**\[experimental\]** Aggregates the per-cell ScType scores into one cell
type call per cluster, see Ianevski et al. (2022).

## Usage

``` r
rs_sc_type_cluster_assignment(sc_type_res, cluster_labels)
```

## Arguments

- sc_type_res:

  List. The ScType results, see
  [`rs_sc_type()`](https://gregorlueg.github.io/bixverse/reference/rs_sc_type.md).

- cluster_labels:

  Integer vector. Cluster assignment per scored cell.

## Value

A list with

- cluster_id - Integer. The cluster id.

- cell_type - Character. The predicted cell type.

- scores - Numeric. The final score for the cluster.

- n_cells - Integer. The number of cells in the cluster.
