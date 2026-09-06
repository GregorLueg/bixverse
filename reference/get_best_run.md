# Get the best run from a stabilised NMF result

Extracts the run with the lowest reconstruction loss and returns it as a
single-run `NmfResult` for downstream methods.

## Usage

``` r
get_best_run(x)

# S3 method for class 'StabilisedNmfResult'
get_best_run(x)
```

## Arguments

- x:

  `StabilisedNmfResult` object.

## Value

An `NmfResult` containing the W/H of the best run.

## Examples

``` r
# the restart with the lowest reconstruction loss
sc <- demo_single_cells()
res <- stabilised_nmf_sc(sc, k = 5L, n_runs = 5L, .verbose = FALSE)
get_best_run(res)
#> NmfResult (single-run HALS NMF)
#>   Source class:     SingleCells
#>   No genes:         30
#>   No cells:         500
#>   No components:    5
#>   Final loss:       5.035e+04
#>   Iterations:       NA
#>   Converged:        TRUE
#>   Preprocessing:    none

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
