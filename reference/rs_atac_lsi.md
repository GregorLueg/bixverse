# LSI of a TF-IDF store

**\[experimental\]**

## Usage

``` r
rs_atac_lsi(
  f_path_cells,
  f_path_genes,
  cells,
  features,
  n_comp,
  solver,
  scale_embeddings,
  scale_max,
  seed,
  verbose
)
```

## Arguments

- f_path_cells:

  String. Cell-based file.

- f_path_genes:

  String. Gene-based file.

- cells:

  Integer vector. Cells (0-based), output row order.

- features:

  Integer vector. Features to use (0-based).

- n_comp:

  Integer. Number of components.

- solver:

  String. `"auto"`, `"exact"`, `"randomised"` or `"covariance"`.

- scale_embeddings:

  Boolean. Z-score each component.

- scale_max:

  Optional numeric. Clip of the z-scores.

- seed:

  Integer. Seed.

- verbose:

  Integer. `0L` - quiet; `1L` - normal; `2L` - detailed.

## Value

A list with `embedding` (cells x `n_comp`), `loadings` (features x
`n_comp`), `singular_values` and `depth_cor`.
