# Get the ADT sample info

Returns sample information in terms of ADT capture, number of features,
etc.

## Usage

``` r
get_adt_sample_info(x)

# S3 method for class 'ADTCounts'
get_adt_sample_info(x)
```

## Arguments

- x:

  The object from which to get the ADT sample information

## Value

A data.table with the ADT sample information

## Examples

``` r
# per-cell ADT capture: non-zero proteins and library size
adt <- generate_single_cell_test_data_adt()
cell_info <- stats::setNames(
  seq_len(nrow(adt$counts)),
  rownames(adt$counts)
)
adt_clr <- new_adt_counts_clr(adt$counts, cell_info = cell_info)
head(get_adt_sample_info(adt_clr))
#>    cell_idx adt_nnz adt_lib_size
#>       <int>   <num>        <num>
#> 1:        1      14         3296
#> 2:        2      15         2101
#> 3:        3      15         3409
#> 4:        4      14         2090
#> 5:        5      15         3664
#> 6:        6      13         3080
```
