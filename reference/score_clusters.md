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

## Examples

``` r
# aggregate the per cell ScType scores onto a clustering
sc <- demo_single_cells()
markers <- data.table::data.table(
  cell_type = rep(sprintf("cell_type_%i", 1:3), each = 10),
  gene_id = sprintf("gene_%02d", 1:30)
)
cell_markers <- prepare_cell_markers(obj = sc, marker_df = markers)
scores <- calc_sc_type_scores(
  sc,
  cell_marker_list = cell_markers,
  .verbose = FALSE
)
clusters <- fast_cluster_sc(
  sc,
  resolutions = 1.0,
  n_centroids = 30L,
  .verbose = FALSE
)
head(score_clusters(scores, cluster_labels = get_data(clusters)$res_1))
#>    cluster_id   cell_type   scores n_cells
#>         <int>      <char>    <num>   <int>
#> 1:          0 cell_type_2 456.6014     173
#> 2:          1 cell_type_3 464.2848     166
#> 3:          2 cell_type_1 449.7363     161

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
