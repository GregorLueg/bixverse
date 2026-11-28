# Write a single cell experiment to h5ad or AnnData Zarr

Streams the counts out of the cell-based binary file, cell batch by cell
batch, into a spec-compliant h5ad file (or an AnnData Zarr v3 store) via
`bixverse-sc-export`, together with the obs and var tables, dense
embeddings and sparse cell x cell graphs supplied from R. Counts are
stored as `float32`. All inputs are checked before the output is
created.

## Usage

``` r
rs_save_h5ad(
  f_path_cells,
  h5_path,
  cell_indices,
  norm,
  obs_index,
  obs,
  var_index,
  var,
  obsm,
  varm,
  obsp,
  uns_json,
  chunk_size,
  zarr,
  gene_indices
)
```

## Arguments

- f_path_cells:

  String. Path to the `counts_cells.bin` file.

- h5_path:

  String. Path of the h5ad file or Zarr directory to create.

- cell_indices:

  Integer vector. The cells to write (0-indexed!), in the order they
  shall appear in the file.

- norm:

  Boolean. Write the normalised instead of the raw counts.

- obs_index:

  Character vector. One name per cell.

- obs:

  Named list. The obs columns, see the R wrapper for the types.

- var_index:

  Character vector. One name per written gene.

- var:

  Named list. The var columns.

- obsm:

  Named list of numeric matrices with one row per cell.

- varm:

  Named list of numeric matrices with one row per written gene.

- obsp:

  Named list of CSR matrices, each a list with `indptr`, `indices`
  (0-indexed) and `data`.

- uns_json:

  String. JSON object written to `uns`.

- chunk_size:

  Integer. Number of cells per streaming batch.

- zarr:

  Boolean. Write an AnnData Zarr store instead of an h5ad file.

- gene_indices:

  Optional integer vector. The genes to write (0-indexed!), in output
  order. `NULL` writes every gene.

## Value

Invisible `NULL`.
