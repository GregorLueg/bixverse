# Run a function over each group with optional progress reporting

Run a function over each group with optional progress reporting

## Usage

``` r
.run_per_group(groups, per_group_fn, .verbose, label)
```

## Arguments

- groups:

  Named list of cell index vectors.

- per_group_fn:

  Function with signature `(cells, name, inner_verbose)`.

- .verbose:

  Verbosity level passed in from the calling function. One level is
  demoted before being forwarded to `per_group_fn`.

- label:

  Character. Label shown on the progress bar.

## Value

Named list of results, one element per group.
