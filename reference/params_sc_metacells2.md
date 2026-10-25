# Wrapper function for the MetaCells2 per-pile parameters

Wrapper function for the MetaCells2 per-pile parameters

## Usage

``` r
params_sc_metacells2(
  target_metacell_size = 48L,
  min_metacell_size = 12L,
  target_metacell_umis = 160000,
  similarity_method = c("abs_pearson", "pearson"),
  log_data = TRUE,
  value_regularisation = 0.142857142857143,
  min_genes = 100L,
  min_gene_total = 50L,
  min_gene_top3 = 4L,
  min_gene_relative_variance = 0.1,
  downsample_min_samples = 750L,
  knn_k_override = NULL,
  min_knn_k = 30L,
  min_gene_fold_factor = 3,
  max_cell_fraction = 0.25,
  gap_skip_cells = 1L,
  min_robust_size_factor = 0.5,
  min_convincing_gene_fold_factor = 3
)
```

## Arguments

- target_metacell_size:

  Integer. Target number of cells per meta cell. Defaults to `48L`.

- min_metacell_size:

  Integer. Candidate meta cells below this size are dissolved into
  outliers. Defaults to `12L`.

- target_metacell_umis:

  Numeric. Target total UMIs per meta cell. Drives the number of seeds
  and the size of the kNN graph. Defaults to `160000.0`.

- similarity_method:

  String. Cell-cell similarity on the selected genes. One of
  `c("abs_pearson", "pearson")`. Defaults to `"abs_pearson"`.

- log_data:

  Boolean. Correlate `log2(x + value_regularisation)` instead of `x`.
  Defaults to `TRUE`.

- value_regularisation:

  Numeric. Constant added before the optional log transform. Defaults to
  `0.14285714285714285`.

- min_genes:

  Integer. Floor on the number of selected genes per pile. The relative
  variance threshold is relaxed until it is met. Defaults to `100L`.

- min_gene_total:

  Integer. Minimum total UMIs (downsampled) for a gene to be selected.
  Defaults to `50L`.

- min_gene_top3:

  Integer. Minimum third-highest UMI count of a gene to be selected.
  Defaults to `4L`.

- min_gene_relative_variance:

  Numeric. Minimum windowed relative variance for a gene to be selected.
  Defaults to `0.1`.

- downsample_min_samples:

  Integer. Floor for the per-pile downsampling target. Defaults to
  `750L`.

- knn_k_override:

  Integer or `NULL`. Fixed k for the balanced kNN graph. `NULL` derives
  it from `target_metacell_umis` and the library sizes. Defaults to
  `NULL`.

- min_knn_k:

  Integer. Floor on the derived k. Defaults to `30L`.

- min_gene_fold_factor:

  Numeric. Minimum log2 fold gap within a candidate meta cell for a gene
  to mark cells as deviant. Defaults to `3.0`.

- max_cell_fraction:

  Numeric. Upper bound on the fraction of cells marked as deviants.
  Defaults to `0.25`.

- gap_skip_cells:

  Integer. Cells skipped over when measuring a gap. Defaults to `1L`.

- min_robust_size_factor:

  Numeric. Multiplier on the target size (and UMIs) above which a meta
  cell is kept unconditionally. Defaults to `0.5`.

- min_convincing_gene_fold_factor:

  Numeric. Minimum log2 fold of any gene for a meta cell below the
  robust size to be kept anyway. Defaults to `3.0`.

## Value

A named list with the following elements:

- target_metacell_size - Integer. Target number of cells per meta cell.
  Defaults to `48L`.

- min_metacell_size - Integer. Candidate meta cells below this size are
  dissolved into outliers. Defaults to `12L`.

- target_metacell_umis - Numeric. Target total UMIs per meta cell.
  Drives the number of seeds and the size of the kNN graph. Defaults to
  `160000.0`.

- similarity_method - String. Cell-cell similarity on the selected
  genes. One of `c("abs_pearson", "pearson")`. Defaults to
  `"abs_pearson"`.

- log_data - Boolean. Correlate `log2(x + value_regularisation)` instead
  of `x`. Defaults to `TRUE`.

- value_regularisation - Numeric. Constant added before the optional log
  transform. Defaults to `0.14285714285714285`.

- min_genes - Integer. Floor on the number of selected genes per pile.
  The relative variance threshold is relaxed until it is met. Defaults
  to `100L`.

- min_gene_total - Integer. Minimum total UMIs (downsampled) for a gene
  to be selected. Defaults to `50L`.

- min_gene_top3 - Integer. Minimum third-highest UMI count of a gene to
  be selected. Defaults to `4L`.

- min_gene_relative_variance - Numeric. Minimum windowed relative
  variance for a gene to be selected. Defaults to `0.1`.

- downsample_min_samples - Integer. Floor for the per-pile downsampling
  target. Defaults to `750L`.

- knn_k_override - Integer or `NULL`. Fixed k for the balanced kNN
  graph. `NULL` derives it from `target_metacell_umis` and the library
  sizes. Defaults to `NULL`.

- min_knn_k - Integer. Floor on the derived k. Defaults to `30L`.

- min_gene_fold_factor - Numeric. Minimum log2 fold gap within a
  candidate meta cell for a gene to mark cells as deviant. Defaults to
  `3.0`.

- max_cell_fraction - Numeric. Upper bound on the fraction of cells
  marked as deviants. Defaults to `0.25`.

- gap_skip_cells - Integer. Cells skipped over when measuring a gap.
  Defaults to `1L`.

- min_robust_size_factor - Numeric. Multiplier on the target size (and
  UMIs) above which a meta cell is kept unconditionally. Defaults to
  `0.5`.

- min_convincing_gene_fold_factor - Numeric. Minimum log2 fold of any
  gene for a meta cell below the robust size to be kept anyway. Defaults
  to `3.0`.

## References

Ben-Kiki, et al., Genome Biol., 2022
