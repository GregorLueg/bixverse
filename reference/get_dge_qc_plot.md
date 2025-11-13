# Return QC plots

Getter function to extract the QC plots from the
[`bulk_dge()`](bulk_dge.md) class. These are added when you run for
example [`qc_bulk_dge()`](qc_bulk_dge.md) and
[`normalise_bulk_dge()`](normalise_bulk_dge.md). You can either leave
the plot choice as `NULL` and provide input when prompted, or you
provide the name. The possible plots that might be in the class

- p1_nb_genes_cohort Proportion of non-zero genes for the samples in the
  respective cohorts (added after using
  [`qc_bulk_dge()`](qc_bulk_dge.md)).

- p2_outliers An outlier plot based on the data from p1, added after
  using [`qc_bulk_dge()`](qc_bulk_dge.md).

- p3_voom_normalization Initial Voom normalisation plot after filtering
  lowly expressed genes. Added after using
  [`normalise_bulk_dge()`](normalise_bulk_dge.md).

- p4_boxplot_normalization Expression levels after normalisation. Added
  after using [`normalise_bulk_dge()`](normalise_bulk_dge.md).

- p5_pca_case_control A PCA plot with the chosen case control category.
  Added if [`calculate_pca_bulk_dge()`](calculate_pca_bulk_dge.md) is
  run.

- p6_batch_correction_plot A PCA plot pre and post batch correction with
  the case-control category overlayed. Added if
  [`batch_correction_bulk_dge()`](batch_correction_bulk_dge.md) is run.

## Usage

``` r
get_dge_qc_plot(object, plot_choice = NULL)
```

## Arguments

- object:

  `bulk_dge` class.

- plot_choice:

  Optional string or integer. Index or name of the plate.

## Value

Returns the DGEList stored in the class.
