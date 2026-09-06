# Assign cell types per cell based on ScType

The cluster-level ScType path
([`score_clusters()`](https://gregorlueg.github.io/bixverse/reference/score_clusters.md))
stamps one label on every cell of a cluster, so a minority population
sharing a Leiden community with a bigger one gets absorbed without
anything flagging it. This runs the scoring per cell instead: the score
matrix is smoothed over the sNN graph via label spreading (Zhou et al.),
then each cell takes its own argmax. Cells whose best score falls below
`score_floor` come back as `NA`.

Pass `cluster_col` to also get the per-cluster composition (purity,
entropy, runner-up cell type) and the hybrid assignment, where clusters
at or above `purity_threshold` keep the blanket cluster-level call and
mixed clusters fall back to the per-cell calls.

## Usage

``` r
assign_sc_type(
  object,
  sc_type_res,
  cluster_col = NULL,
  sctype_cell_params = params_sctype_cells(),
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells`.

- sc_type_res:

  `ScTypeResults`, see
  [`calc_sc_type_scores()`](https://gregorlueg.github.io/bixverse/reference/calc_sc_type_scores.md).

- cluster_col:

  Optional string. Name of the obs column with the cluster assignment.
  If provided, the composition and hybrid assignment are returned on top
  of the per-cell calls.

- sctype_cell_params:

  List. Output of
  [`params_sctype_cells()`](https://gregorlueg.github.io/bixverse/reference/params_sctype_cells.md).

- .verbose:

  Boolean or integer. Controls verbosity.

## Value

An `ScTypeCellResults` results class.

## References

Ianevski et al., Nat Comm, 2022. Zhou et al., NIPS, 2004.

## Examples

``` r
# per cell calls from the ScType scores
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
res <- assign_sc_type(sc, sc_type_res = scores, .verbose = FALSE)
table(res$assignments)
#> 
#> cell_type_1 cell_type_2 cell_type_3 
#>         167         166         167 

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
