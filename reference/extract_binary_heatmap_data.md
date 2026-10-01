# Extract plot-ready data for a binary heatmap

Filters, orders and optionally bins a logical samples x features matrix,
e.g. the regulon on/off calls from
[`binarise_regulon_activity()`](https://gregorlueg.github.io/bixverse/reference/binarise_regulon_activity.md),
so it can be drawn as a single raster with
`bixverse.plots::plot_binary_heatmap()`.

Features are clustered within their group on the Jaccard distance, so
shared absences do not count as similarity. Samples are clustered within
their group on the Hamming distance. Groups larger than `max_cluster_n`
are instead ordered by barycentre, the mean plot position of the
features a sample has on, which avoids the quadratic distance matrix.
With more samples than `max_cols`, consecutive samples within a group
are collapsed into roughly `max_cols` bins that hold the fraction of
samples on. Every group keeps at least one bin.

## Usage

``` r
extract_binary_heatmap_data(
  binary_mat,
  sample_groups = NULL,
  feature_groups = NULL,
  heatmap_params = params_binary_heatmap(),
  .verbose = TRUE
)
```

## Arguments

- binary_mat:

  Logical matrix of samples x features with unique row and column names.

- sample_groups:

  Optional named character vector or factor mapping every row name of
  `binary_mat` to a group. Factor levels set the group order, otherwise
  groups are sorted.

- feature_groups:

  Optional named character vector or factor mapping every column name of
  `binary_mat` to a group. Factor levels set the group order, otherwise
  groups are sorted.

- heatmap_params:

  List. Output of
  [`params_binary_heatmap()`](https://gregorlueg.github.io/bixverse/reference/params_binary_heatmap.md).

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

An object of class `BinaryHeatmapData`, a list with:

- mat - Numeric matrix of features x columns in plot order. Values are
  0/1, or the fraction of samples on if binned.

- col_annot - data.table with one row per column: `col_idx`, `group`
  (factor) and `n_samples` in that column.

- row_annot - data.table with one row per feature: `row_idx`, `feature`
  and `group` (factor).

- binned - Boolean. Whether the samples were binned.

Without groups, `group` is a single level `"all"`.

## References

Aibar, et al., Nat Methods, 2017

## Examples

``` r
set.seed(7L)
binary_mat <- matrix(runif(200 * 30) > 0.7, nrow = 200)
dimnames(binary_mat) <- list(
  sprintf("cell_%i", 1:200),
  sprintf("regulon_%i", 1:30)
)
sample_groups <- setNames(rep(c("a", "b"), each = 100), rownames(binary_mat))
res <- extract_binary_heatmap_data(
  binary_mat,
  sample_groups = sample_groups,
  .verbose = FALSE
)
dim(res$mat)
#> [1]  30 200
```
