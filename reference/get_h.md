# Get the H (cell activations) matrix

Get the H (cell activations) matrix

## Usage

``` r
get_h(x)

# S3 method for class 'NmfResult'
get_h(x)

# S3 method for class 'StabilisedNmfResult'
get_h(x)

# S3 method for class 'ConsensusNmfResult'
get_h(x)
```

## Arguments

- x:

  An object holding NMF results.

## Examples

``` r
# per cell activations of a five factor NMF
sc <- demo_single_cells()
res <- nmf_sc(sc, k = 5L, .verbose = FALSE)
dim(get_h(res))
#> [1]   5 500

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
