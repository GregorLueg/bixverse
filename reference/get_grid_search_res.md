# Get the grid search results

Getter function to extract the grid search results. If not found will
return `NULL`.

## Usage

``` r
get_grid_search_res(object)
```

## Arguments

- object:

  The class, see [`bulk_coexp()`](bulk_coexp.md).

## Value

data.table with the grid search results (if found. Otherwise `NULL`.)
