# Get the ready obs data from various sub method

Helper method that creates data.tables with cell indices which were used
in the given analysis + the values that are to be added to the obs table
in the DuckDB.

## Usage

``` r
get_data(x, columns = NULL, ...)

# S3 method for class 'CellQc'
get_data(x, ...)

# S3 method for class 'ScListRes'
get_data(x, ...)

# S3 method for class 'ScMatrixRes'
get_data(x, columns = NULL, ...)

# S3 method for class 'ScrubletRes'
get_data(x, ...)

# S3 method for class 'BoostRes'
get_data(x, ...)

# S3 method for class 'ScDblFinderRes'
get_data(x, ...)

# S3 method for class 'SingleCellFastClusters'
get_data(x, ...)

# S3 method for class 'NmfResult'
get_data(x, ...)

# S3 method for class 'ConsensusNmfResult'
get_data(x, ...)
```

## Arguments

- x:

  An object to get the data from.

- columns:

  Optional string. For some of the functions you can decide to only
  extract specific columns.

- ...:

  Other parameters

## Value

Returns a data.table with a cell_idx column for the cells included in
the analysis and additional columns to be added to the obs table.

## Examples

``` r
# the cell indices and cluster memberships a run produced
sc <- demo_single_cells()
res <- fast_cluster_sc(sc, resolutions = 1.0, .verbose = FALSE)
head(get_data(res), 3)
#>    cell_idx res_1
#>       <int> <int>
#> 1:        1     2
#> 2:        2     0
#> 3:        3     1

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
