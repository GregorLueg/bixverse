# Seurat label transfer

Seurat's `FindTransferAnchors(reduction = "cca")` between the kept cells
of both experiments on shared genes, then `TransferData()`: the anchor
labels weighted per query cell in `weight_embedding`. For RNA to ATAC
the query is the gene activity from
[`gene_activity_sc()`](https://gregorlueg.github.io/bixverse/reference/gene_activity_sc.md)
and the weights come from the ATAC LSI, as in the Signac vignette; RNA
to RNA works the same with the query PCA.

The predictions and their confidence are written to the query's obs as
`predicted_<label_column>` and `confidence_<label_column>`.

## Usage

``` r
transfer_labels_seurat_sc(
  reference,
  query,
  label_column,
  weight_embedding,
  features = NULL,
  transfer_params = params_sc_seurat_transfer(),
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- reference:

  `SingleCells` class. The annotated experiment (RNA).

- query:

  `SingleCells` class. The experiment to label (RNA or gene activity),
  matched to the reference on `gene_id`.

- label_column:

  String. Reference obs column with the labels. Missing labels form a
  level of their own.

- weight_embedding:

  Numeric matrix. Kept query cells x dimensions for the weights, e.g.
  the `"lsi"` embedding of the ATAC side. Rows are matched to the query
  by position; with row names, those are checked against the kept query
  cells.

- features:

  Optional character vector. Gene ids to anchor on. `NULL` uses the
  reference HVGs. Only those in both experiments are used.

- transfer_params:

  List. Output of
  [`params_sc_seurat_transfer()`](https://gregorlueg.github.io/bixverse/reference/params_sc_seurat_transfer.md).

- seed:

  Integer. Seed of the kNN backend.

- .verbose:

  Boolean or integer. Controls verbosity.

## Value

A data.table, one row per kept query cell, with `cell_id`,
`predicted_<label_column>`, `confidence_<label_column>` (the top score)
and one `score_<label>` column per label.

## References

Stuart, et al., Cell, 2019
