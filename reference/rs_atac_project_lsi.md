# Project cells of a TF-IDF store onto a fitted LSI

**\[experimental\]** Seurat's `ProjectSVD(mode = "lsi")`: the query's
own TF-IDF times the loadings over the singular values, no centring or
scaling.

## Usage

``` r
rs_atac_project_lsi(
  f_path_cells,
  f_path_genes,
  cells,
  features,
  loadings,
  singular_values
)
```

## Arguments

- f_path_cells:

  String. Cell-based file of the query.

- f_path_genes:

  String. Gene-based file of the query (TF-IDF norm).

- cells:

  Integer vector. Query cells (0-based), output row order.

- features:

  Integer vector. Query features (0-based), in the row order of
  `loadings`.

- loadings:

  Numeric matrix. Features x components of the fit.

- singular_values:

  Numeric vector. Singular values of the fit.

## Value

Numeric matrix, cells x components.
