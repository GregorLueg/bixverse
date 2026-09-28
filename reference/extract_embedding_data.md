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

## Examples

``` r
# PCA coordinates with a cell annotation riding along
sc <- demo_single_cells()
dt <- extract_embedding_data(sc, "pca", obs_cols = "cell_grp")
head(dt[, c("cell_id", "dim_1", "dim_2", "cell_grp")])
#>     cell_id      dim_1      dim_2    cell_grp
#>      <char>      <num>      <num>      <char>
#> 1: cell_001 -0.3694211  3.0682013 cell_type_1
#> 2: cell_002  2.2331283 -2.5044506 cell_type_2
#> 3: cell_003 -2.4834063  0.5060492 cell_type_3
#> 4: cell_004  1.3577558  2.6232290 cell_type_1
#> 5: cell_005  2.1342576 -0.2906672 cell_type_2
#> 6: cell_006 -2.6232185  0.5387377 cell_type_3

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
