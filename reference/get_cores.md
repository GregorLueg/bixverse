# Helper function for core detection.

Identifies the number of cores/workers to use for parallel tasks. Will
default to half of the available cores, to a maximum of
`abs_max_workers`.

## Usage

``` r
get_cores(abs_max_workers = 8L)
```

## Arguments

- abs_max_workers:

  Integer. Absolute maximum number of cores to use. Defaults to `8L`.

## Value

Integer. Number of cores to use.
