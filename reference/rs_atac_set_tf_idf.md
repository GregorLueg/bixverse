# Switch a store to TF-IDF over a set of cells

**\[experimental\]** Header-only change of the norm; the chunks are not
rewritten.

## Usage

``` r
rs_atac_set_tf_idf(f_path_cells, f_path_genes, cells, scale)
```

## Arguments

- f_path_cells:

  String. Cell-based file.

- f_path_genes:

  String. Gene-based file.

- cells:

  Optional integer vector. Cells for the IDF (0-based), `NULL` for all.

- scale:

  Numeric. Scale factor.

## Value

The IDF per feature.
