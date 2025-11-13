# Helper function to prepare data for GSEA

Helper function to prepare data for GSEA

## Usage

``` r
prep_stats_pathways(stats, pathways, min_size = 5L, max_size = 500L)
```

## Arguments

- stats:

  Named numeric vector. The gene level statistic.

- pathways:

  List. A named list with each element containing the genes for this
  pathway.

- min_size:

  Integer. Minimum size of the gene sets. Defaults to `5L`.

- max_size:

  Integer. Maximum size of the gene sets. Defaults to `500L`.

## Value

A list with the following elements:

- stats - The sorted gene level statistics

- pathways - The prepared pathways for further usage in the functions

- pathway_sizes - The sizes of the pathways
