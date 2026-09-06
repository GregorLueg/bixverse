# Print a BoostRes object

Print a BoostRes object

## Usage

``` r
# S3 method for class 'BoostRes'
print(x, ...)
```

## Arguments

- x:

  A `BoostRes` object.

- ...:

  Ignored.

## Value

Invisible `x`.

## Examples

``` r
# doublet calls from the boosted classifier
sc <- demo_single_cells(prepped = FALSE)
res <- doublet_detection_boost_sc(
  sc,
  boost_params = params_boost(
    hvg = list(min_gene_var_pctl = 0.0),
    pca = list(no_pcs = 10L),
    n_iters = 5L
  ),
  .verbose = FALSE
)
print(res)
#> BoostRes: 500 cells, 1 doublets (0.2%)
#>   Score range: [0.0125, 0.8289]

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
