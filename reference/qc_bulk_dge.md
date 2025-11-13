# QC on the bulk dge data

This function will do QC on the bulk data and remove outlier samples
that show substantially lower expression of genes compared to the rest
of the data and remove lowly expressed genes.

## Usage

``` r
qc_bulk_dge(
  object,
  group_col,
  outlier_threshold = 2,
  min_prop = 0.2,
  min_count = 10,
  .verbose = TRUE
)
```

## Arguments

- object:

  The underlying class, see [`bulk_dge()`](bulk_dge.md).

- group_col:

  String. The column in the metadata that will contain the contrast
  groups. Needs to be part of the metadata stored in the class.

- outlier_threshold:

  Float. Number of standard deviations in terms of percentage genes
  detected you allow before removing a sample. Defaults to `2`.

- min_prop:

  Float. Minimum proportion of samples in which the gene has to be
  identified in.

- min_count:

  Float. Minimum number of counts (cpm) to be detected in min_prop of
  the samples (in cohorts defined by groups_coll)

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

Returns the class with the `processed_data` data slot populated and
applied parameters added to the `params` slot.
