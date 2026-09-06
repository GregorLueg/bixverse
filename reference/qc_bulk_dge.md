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

  The underlying class, see
  [`BulkDge()`](https://gregorlueg.github.io/bixverse/reference/BulkDge.md).

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

## Examples

``` r
# drop outlier samples and lowly expressed genes
syn <- synthetic_bulk_cor_matrix()
meta <- data.table::data.table(
  sample_id = colnames(syn$counts),
  case_control = rep(c("case", "control"), each = 50)
)
object <- BulkDge(raw_counts = syn$counts, meta_data = meta)
object <- qc_bulk_dge(object, group_col = "case_control", .verbose = FALSE)
dim(get_outputs(object)$raw_counts_filtered)
#> [1] 1000   98
```
