# Rescue MAD-flagged cells that fall within safe bounds

For each metric in `rescue_thresholds`, cells whose value is within
`[lower, upper]` are un-flagged from MAD. Hard-threshold flags are not
cleared; a warning is issued when overlaps occur. Calling `rescue_cells`
overwrites any prior rescue set for the metrics named.

## Usage

``` r
rescue_cells(x, rescue_thresholds)
```

## Arguments

- x:

  A `CellQc` object.

- rescue_thresholds:

  Named list of numeric vectors with `lower` and/or `upper`.

## Value

Updated `CellQc`.
