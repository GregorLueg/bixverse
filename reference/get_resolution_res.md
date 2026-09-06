# Return the resolution results

Getter function to get the resolution results (if available).

## Usage

``` r
get_resolution_res(object)
```

## Arguments

- object:

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

## Value

If resolution results were found, returns the data.table. Otherwise,
throws a warning and returns NULL.

## Examples

``` r
# cluster counts and modularity across the Leiden resolution sweep
syn <- synthetic_bulk_cor_matrix()
mat <- log1p(t(syn$counts))
meta <- data.table::data.table(sample_id = rownames(mat))
object <- BulkCoExp(raw_data = mat, meta_data = meta)
object <- preprocess_bulk_coexp(object, hvg = 200L, .verbose = FALSE)
object <- cor_module_processing(
  object,
  cor_method = "pearson",
  .verbose = FALSE
)
object <- cor_module_graph_check_res(
  object,
  parallel = FALSE,
  .verbose = FALSE
)
head(get_resolution_res(object))
#> Key: <resolution>
#>    resolution no_clusters modularity good_clusters avg_size max_size
#>         <num>       <int>      <num>         <int>    <num>    <int>
#> 1:  0.1000000           3  0.6603861             3 61.33333       65
#> 2:  0.1389495           3  0.6603861             3 61.33333       65
#> 3:  0.1930698           3  0.6603861             3 61.33333       65
#> 4:  0.2682696           3  0.6603861             3 61.33333       65
#> 5:  0.3727594           3  0.6603861             3 61.33333       65
#> 6:  0.5179475           3  0.6603861             3 61.33333       65
```
