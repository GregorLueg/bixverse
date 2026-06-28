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
