# Print a ScDblFinderRes object

Print a ScDblFinderRes object

## Usage

``` r
# S3 method for class 'ScDblFinderRes'
print(x, ...)
```

## Arguments

- x:

  A `ScDblFinderRes` object.

- ...:

  Ignored.

## Value

Invisible `x`.

## Examples

``` r
# threshold, score range and the cluster count behind the calls
sc <- demo_single_cells(prepped = FALSE)
res <- scdblfinder_sc(
  sc,
  scdblfinder_params = params_scdblfinder(
    pca = list(no_pcs = 10L),
    n_genes = 25L,
    cxds_genes = 25L
  ),
  .verbose = FALSE
)
print(res)
#> ScDblFinderRes: 500 cells, 16 doublets (3.2%)
#>   Threshold:        0.4517
#>   Score range:      [0.0361, 0.9481]
#>   Final clusters:   3
#>   Features available: FALSE

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
