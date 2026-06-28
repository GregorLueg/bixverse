# Update doublet calls and summary statistics for a new Scrublet threshold

Recomputes `predicted_doublets`, `z_scores`, and the three rate fields
(`detected_doublet_rate`, `detectable_doublet_fraction`,
`overall_doublet_rate`) for either the full result or a single group.

## Usage

``` r
.update_scrublet_threshold(scrublet_res, threshold, sample_name, .verbose)
```

## Arguments

- scrublet_res:

  A `ScrubletRes` object.

- threshold:

  Numeric. The new score threshold to apply.

- sample_name:

  Character or `NULL`. If `NULL`, updates apply to the whole object
  (ungrouped). Otherwise, updates are scoped to the named group.

- .verbose:

  Logical. If `TRUE`, prints updated rate summaries to the console.

## Value

The `ScrubletRes` object with updated doublet calls and summary
statistics.
