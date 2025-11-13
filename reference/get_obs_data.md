# Get the ready obs data from various sub method

Helper method that creates data.tables with cell indices which were used
in the given analysis + the values that are to be added to the obs table
in the DuckDB.

## Usage

``` r
get_obs_data(x, ...)

# S3 method for class 'sc_proportion_res'
get_obs_data(x, ...)

# S3 method for class 'scrublet_res'
get_obs_data(x, ...)

# S3 method for class 'boost_res'
get_obs_data(x, ...)
```

## Arguments

- x:

  An object to set gene mapping for

- ...:

  Other parameters

## Value

Returns a data.table with a cell_idx column for the cells included in
the analysis and additional columns to be added to the obs table.
