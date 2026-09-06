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

## Examples

``` r
# cluster memberships carry their own cell_idx, so they join straight on
sc <- demo_single_cells()
clusters <- fast_cluster_sc(sc, resolutions = 1.0, .verbose = FALSE)
sc <- add_sc_new_obs(sc, get_data(clusters))
head(get_sc_obs(sc), 3)
#>    cell_idx  cell_id    cell_grp batch_index   nnz lib_size to_keep res_1
#>       <int>   <char>      <char>       <num> <num>    <num>  <lgcl> <int>
#> 1:        1 cell_001 cell_type_1           1    34      278    TRUE     2
#> 2:        2 cell_002 cell_type_2           1    31      333    TRUE     0
#> 3:        3 cell_003 cell_type_3           1    32      413    TRUE     1

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
