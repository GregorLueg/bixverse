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
  weights if `sensitivity = TRUE`. If not provided, defaults to `0` as
  in the original implementation.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

An `ScTypeResults` results class
