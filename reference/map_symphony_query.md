# Map a SingleCells query onto a Symphony reference

Projects query cells through the reference's PCA loadings and applies
the cached MoE batch correction. Gene matching is by name against the
reference's stored HVG names; unmatched HVGs become zero columns.
Results are attached to the query as embeddings.

## Usage

``` r
map_symphony_query(
  reference,
  query,
  batch_column = NULL,
  additional_batch_columns = NULL,
  params = params_symphony_map(),
  .verbose = TRUE
)
```

## Arguments

- reference:

  `SymphonyReference`.

- query:

  `SingleCells` query.

- batch_column:

  Optional string. If `NULL`, no batch correction is applied (z_corr =
  z_pca).

- additional_batch_columns:

  Optional character vector.

- params:

  List. Output of
  [`params_symphony_map()`](https://gregorlueg.github.io/bixverse/reference/params_symphony_map.md).

- .verbose:

  Boolean or integer.

## Value

The `query` object with embeddings `"symphony"` (z_corr),
`"symphony_pca"` (z_pca) and `"symphony_r"` (soft cluster assignments,
transposed to N_q x K).

## Examples

``` r
# project a second data set into the reference embedding
ref <- demo_single_cells(
  syn_data_params = params_sc_synthetic_data(
    n_cells = 500L,
    n_genes = 50L,
    n_batches = 2L
  )
)
symphony_ref <- build_symphony_ref(
  ref,
  batch_column = "batch_index",
  hvg = get_hvg(ref) + 1L,
  harmony_params = params_sc_harmony(k = 10L),
  no_pcs = 10L,
  .verbose = FALSE
)
query <- demo_single_cells(prepped = FALSE, seed = 7L)
query <- map_symphony_query(symphony_ref, query = query, .verbose = FALSE)
dim(get_embedding(query, "symphony"))
#> [1] 500  10

unlink(c(ref@dir_data, query@dir_data), recursive = TRUE, force = TRUE)
```
