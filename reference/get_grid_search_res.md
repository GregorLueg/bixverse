# Get the grid search results

Getter function to extract the grid search results. If not found will
return `NULL`.

## Usage

``` r
get_grid_search_res(object)
```

## Arguments

- object:

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

## Value

data.table with the grid search results (if found. Otherwise `NULL`.)

## Examples

``` r
# the grid search table
syn <- generate_gene_module_data(n_samples = 24L, n_genes = 60L)
obj <- BulkCoExp(syn$data, syn$meta_data)
obj <- preprocess_bulk_coexp(obj, hvg = NULL, .verbose = FALSE)
obj <- dgrdl_grid_search(
  obj,
  neighbours_vec = c(3L, 5L),
  dict_size_vec = c(4L, 6L),
  seed_vec = 123L,
  .verbose = FALSE
)
get_grid_search_res(obj)
#>     seed dict_size k_neighbours reconstruction_errs feature_laplacian_objective
#>    <num>     <num>        <num>               <num>                       <num>
#> 1:   123         4            3            7.809673                  0.08282019
#> 2:   123         4            5            8.081732                  0.06261729
#> 3:   123         6            3            3.477399                  0.09445959
#> 4:   123         6            5            3.794702                  0.07474547
#>    sample_laplacian_objective
#>                         <num>
#> 1:                  0.3560355
#> 2:                  1.7663868
#> 3:                  0.5568773
#> 4:                  2.7064393
```
