# Helper function to manually readjust Scrublet thresholds

Can update the Scrublet thresholding after manual introspection of the
histogram.

## Usage

``` r
call_doublets_manual(scrublet_res, threshold, .verbose = TRUE)
```

## Arguments

- scrublet_res:

  `scrublet_res` result class.

- threshold:

  Numeric. The new threshold to use.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

`scrublet_res` class with updated doublet calls based on the new
threshold.
