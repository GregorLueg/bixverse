# Get scores

Get scores

## Usage

``` r
get_scores(x, ...)

# S3 method for class 'ScDblFinderRes'
get_scores(x, ..., score_type = c("weighted", "cxds_scores"))

# S3 method for class 'ScTypeResults'
get_scores(x, ...)

# S3 method for class 'ScTypeCellResults'
get_scores(x, ...)
```

## Arguments

- x:

  An object to get scores from.

- ...:

  Additional arguments passed to methods.

- score_type:

  Either `"weighted"` or `"cxds_scores"`.

## Value

The score matrix or data.table held by the object, depending on the
class. Methods exist for `ScTypeResults`, `ScTypeCellResults` and
`ScDblFinderRes`.

A numeric matrix of cells x cell types.

A numeric vector with the winning score per cell.

## Examples

``` r
# the ScType score matrix of cell types by cluster
sc <- demo_single_cells()
markers <- data.table::data.table(
  cell_type = rep(sprintf("cell_type_%i", 1:3), each = 10),
  gene_id = sprintf("gene_%02d", 1:30)
)
cell_markers <- prepare_cell_markers(obj = sc, marker_df = markers)
res <- calc_sc_type_scores(sc, cell_marker_list = cell_markers,
                           .verbose = FALSE)
dim(get_scores(res))
#> [1] 500   3

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
