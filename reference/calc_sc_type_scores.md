# Calculate ScType scores per cell

Implements the approach from

## Usage

``` r
calc_sc_type_scores(
  object,
  cell_marker_list,
  sensitivity = TRUE,
  weight_floor = NULL,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` or `SingleCellsMultiModal`

- cell_marker_list:

  A list, see
  [`prepare_cell_markers()`](https://gregorlueg.github.io/bixverse/reference/prepare_cell_markers.md).

- sensitivity:

  Boolean. Shall shared marker genes be downweighted (like in the
  original reference). Defaults to `TRUE`.

- weight_floor:

  Optional numeric. A value between 0 to 1 and sets the floor for the
  weights if `sensitivity = TRUE`. If not provided, defaults to `0.1`.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

An `ScTypeResults` results class

## Examples

``` r
# ScType scores against the planted marker genes
sc <- demo_single_cells()
markers <- data.table::data.table(
  cell_type = rep(sprintf("cell_type_%i", 1:3), each = 10),
  gene_id = sprintf("gene_%02d", 1:30)
)
cell_markers <- prepare_cell_markers(obj = sc, marker_df = markers)
res <- calc_sc_type_scores(
  sc,
  cell_marker_list = cell_markers,
  .verbose = FALSE
)
dim(get_scores(res))
#> [1] 500   3

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
