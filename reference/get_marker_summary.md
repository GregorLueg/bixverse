# Get the per-gene marker summaries across all rivals

Returns the summaries a marker is judged on: the AUROC of the reference
against its rivals reduced to one row per gene and reference group. Rank
on `median_auroc` for a marker that survives a single closely related
rival, or on `min_auroc` when it has to beat every rival unambiguously.

## Usage

``` r
get_marker_summary(x)

# S3 method for class 'ScSpecificMarkers'
get_marker_summary(x)
```

## Arguments

- x:

  `ScSpecificMarkers` object.

## Value

A copy of the summary data.table.
