# Returns the available features for single cell applications

Returns a data.table with available features in the obs table and in the
count matrices.

## Usage

``` r
get_sc_available_features(object)
```

## Arguments

- object:

  `SingleCells`, `MetaCells` (or potentially other) class.

## Value

A data.table with available features.

## Examples

``` r
# what can be queried from the obs table and the counts
sc <- demo_single_cells(prepped = FALSE)
head(get_sc_available_features(sc), 3)
#>    feature_name origin
#>          <char> <char>
#> 1:     cell_idx    obs
#> 2:      cell_id    obs
#> 3:     cell_grp    obs

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
