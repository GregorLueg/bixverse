# Add an obs table derived from a method to the SingleCells.

Add an obs table derived from a method to the SingleCells.

## Usage

``` r
add_sc_new_obs(object, obs_data)
```

## Arguments

- object:

  `SingleCells` class.

- obs_data:

  data.table. A data.table you generated with
  [`get_data()`](https://gregorlueg.github.io/bixverse/reference/get_data.md)
  on some sub class.

## Value

The class with updated obs table in the DuckDB
