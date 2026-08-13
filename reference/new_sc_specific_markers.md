# Constructor for the one-vs-many specific marker results

Holds the results of
[`find_specific_markers_sc()`](https://gregorlueg.github.io/bixverse/reference/find_specific_markers_sc.md):
the statistics of each reference group against each of its rivals
separately, plus the per-gene summaries across those rivals. The
summaries are the interesting part, as a gene only marks the reference
if it holds up against every rival.

## Usage

``` r
new_sc_specific_markers(summary, per_comparison, params)
```

## Arguments

- summary:

  data.table. The per-gene summaries across all rivals.

- per_comparison:

  data.table. The per-rival statistics.

- params:

  List. The parameters the run used.

## Value

Generates the `ScSpecificMarkers` class.
