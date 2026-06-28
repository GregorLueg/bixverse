# Get scores

Get scores

## Usage

``` r
get_scores(x, ...)

# S3 method for class 'ScDblFinderRes'
get_scores(x, ..., score_type = c("weighted", "cxds_scores"))

# S3 method for class 'ScTypeResults'
get_scores(x, ...)
```

## Arguments

- x:

  An object to get scores from.

- ...:

  Additional arguments passed to methods.

- score_type:

  Either `"weighted"` or `"cxds_scores"`.

## Value

A numeric matrix of cells x cell types.
