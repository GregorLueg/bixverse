# Extract embedding coordinates for plotting

Pulls an embedding into a long data.table with standardised coordinate
columns (`dim_1`, `dim_2`, ...) and, optionally, observation metadata.
The embedding name is stored as an `embedding` attribute for axis
labelling.

## Usage

``` r
extract_embedding_data(object, embedding, obs_cols = NULL, ...)
```

## Arguments

- object:

  A single cell class.

- embedding:

  String. Name of the embedding (e.g. `"umap"`, `"pca"`).

- obs_cols:

  Optional character vector. Obs columns to attach.

- ...:

  Additional arguments forwarded to
  [`get_embedding()`](https://gregorlueg.github.io/bixverse/reference/get_embedding.md)
  (e.g. `modality`).

## Value

A data.table with `cell_id`, `dim_*` columns and any requested obs
columns.
