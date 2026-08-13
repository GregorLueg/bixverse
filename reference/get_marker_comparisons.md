# Get the per-rival marker statistics

Returns the statistics that the summaries of
[`get_marker_summary()`](https://gregorlueg.github.io/bixverse/reference/get_marker_summary.md)
are built from, i.e., one row per gene, reference group and rival.
Useful to find out which rival a gene fails against.

## Usage

``` r
get_marker_comparisons(x)

# S3 method for class 'ScSpecificMarkers'
get_marker_comparisons(x)
```

## Arguments

- x:

  `ScSpecificMarkers` object.

## Value

A copy of the per comparison data.table.
