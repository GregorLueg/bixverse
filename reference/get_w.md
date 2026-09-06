# Get the W (gene loadings) matrix

Get the W (gene loadings) matrix

## Usage

``` r
get_w(x)

# S3 method for class 'NmfResult'
get_w(x)

# S3 method for class 'StabilisedNmfResult'
get_w(x)

# S3 method for class 'ConsensusNmfResult'
get_w(x)
```

## Arguments

- x:

  An object holding NMF results.

## Examples

``` r
# gene loadings of a five factor NMF
sc <- demo_single_cells()
res <- nmf_sc(sc, k = 5L, .verbose = FALSE)
dim(get_w(res))
#> [1] 30  5

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
