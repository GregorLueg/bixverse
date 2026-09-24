# Calculate cell type average silhouette width

Average silhouette width on cell type labels in the embedding, rescaled
to `[0, 1]` via `(s + 1) / 2` as in scIB. Higher values mean cell types
stay separated after correction. Counterpart to
[`calculate_batch_asw_sc()`](https://gregorlueg.github.io/bixverse/reference/calculate_batch_asw_sc.md)
on the biology side.

## Usage

``` r
calculate_cell_type_asw_sc(
  object,
  cell_type_column,
  embd_to_use = "pca",
  max_cells = 5000L,
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` or `SingleCellsSubset` class.

- cell_type_column:

  String. The column with the cell type labels in the obs data of the
  class.

- embd_to_use:

  String. Which embedding to compute the ASW on. Defaults to `"pca"`.

- max_cells:

  Integer or `NULL`. If not `NULL`, subsample to this many cells for
  performance. Defaults to `5000L`.

- seed:

  Integer. Seed for subsampling reproducibility.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

A `CellTypeAswScores` object with the following elements

- per_cell - Per-cell rescaled silhouette scores in `[0, 1]`.

- mean_asw - Mean rescaled silhouette width.

- median_asw - Median rescaled silhouette width.

- n_cell_types - Number of cell types.

- embedding_used - Which embedding the ASW was computed on.

## References

Luecken, et al., Nat. Methods, 2022

## Examples

``` r
# cell type silhouette width on the PCA embedding
sc <- demo_single_cells(
  syn_data_params = params_sc_synthetic_data(
    n_cells = 600L, n_genes = 50L, n_batches = 3L
  )
)
calculate_cell_type_asw_sc(
  sc,
  cell_type_column = "cell_grp",
  .verbose = FALSE
)
#> Cell Type Silhouette Width (rescaled)
#>   Cells: 600 | Cell types: 3 | Embedding: pca
#>   Mean ASW:    0.5309 (0.5 = no structure, 1 = separated)
#>   Median ASW:  0.5313

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
