# Calculate a summary of integration metrics

Runs the batch mixing and (if cell type labels are given) the biological
conservation metrics in one go and returns one row per call, so results
across correction methods can be `rbind`-ed into one table. Every column
is on `[0, 1]` (PCR comparison can go negative) and higher is better,
the scIB convention:

Batch mixing:

- kbet_accept - `1 - ` kBET rejection rate.

- batch_asw - `mean(1 - |s|)` over the per-cell batch silhouettes.

- ilisi - Normalised iLISI.

- pcr_comparison - `(pre - post) / pre` of the batch PCR.

Biological conservation:

- clisi - Normalised cLISI.

- cell_type_asw - Rescaled cell type silhouette width.

- graph_connectivity - Mean graph connectivity over cell types.

The kNN metrics read the kNN graph currently stored in the object, so
recompute the neighbours on the corrected embedding first. Embedding
metrics are `NA` if `embd_to_use = NULL` (e.g. BBKNN, which only returns
a graph).

## Usage

``` r
calculate_integration_metrics_sc(
  object,
  batch_column,
  cell_type_column = NULL,
  embd_to_use = "pca",
  max_cells = 5000L,
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` or `SingleCellsSubset` class.

- batch_column:

  String. The column with the batch information in the obs data of the
  class.

- cell_type_column:

  Optional string. The column with the cell type labels. If `NULL`, the
  conservation metrics are `NA`.

- embd_to_use:

  Optional string. The embedding for ASW and PCR. Defaults to `"pca"`.

- max_cells:

  Integer or `NULL`. Subsampling for the silhouette widths. Defaults to
  `5000L`.

- seed:

  Integer. Seed for subsampling reproducibility.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

A one-row data.table with the columns `embedding`, `kbet_accept`,
`batch_asw`, `ilisi`, `pcr_comparison`, `clisi`, `cell_type_asw` and
`graph_connectivity`.

## References

Luecken, et al., Nat. Methods, 2022

## Examples

``` r
# all metrics on the uncorrected PCA
sc <- demo_single_cells(
  syn_data_params = params_sc_synthetic_data(
    n_cells = 600L, n_genes = 50L, n_batches = 3L
  )
)
calculate_integration_metrics_sc(
  sc,
  batch_column = "batch_index",
  cell_type_column = "cell_grp",
  .verbose = FALSE
)
#>    embedding kbet_accept batch_asw ilisi pcr_comparison     clisi cell_type_asw
#>       <char>       <num>     <num> <num>          <num>     <num>         <num>
#> 1:       pca   0.2466667 0.9330131   0.4             NA 0.5044248     0.5309407
#>    graph_connectivity
#>                 <num>
#> 1:                  1

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
