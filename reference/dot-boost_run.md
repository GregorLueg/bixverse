# Run boosted doublet detection on a set of cells

Run boosted doublet detection on a set of cells

## Usage

``` r
.boost_run(object, cells_to_use, boost_params, seed, streaming, .verbose)
```

## Arguments

- object:

  A `SingleCells` object.

- cells_to_use:

  Integer vector of 0-indexed cell indices.

- boost_params:

  List of boost parameters from
  [`params_boost()`](https://gregorlueg.github.io/bixverse/reference/params_boost.md).

- seed:

  Integer. Random seed.

- streaming:

  Optional logical. Whether to stream the count data. If `NULL`,
  resolved automatically via
  [`auto_streaming()`](https://gregorlueg.github.io/bixverse/reference/auto_streaming.md).

- .verbose:

  Logical or integer. Verbosity level.

## Value

A `BoostRes` object with `cell_indices` set as an attribute.
