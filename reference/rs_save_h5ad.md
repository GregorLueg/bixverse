# Write a single cell experiment to h5ad

Streams the counts out of the cell-based binary file, cell batch by cell
batch, into a spec-compliant h5ad file via
[scx-core](https://github.com/btraven00/scx), together with the obs and
var tables, dense embeddings and sparse cell x cell graphs supplied from
R. Counts are stored as `float32`. All inputs are checked before the
file is created.

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
  chunk_size
)
```

## Arguments

- f_path_cells:

  String. Path to the `counts_cells.bin` file.

- h5_path:

  String. Path of the h5ad file to create.

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

  Character vector. One name per gene.

- var:

  Named list. The var columns.

- obsm:

  Named list of numeric matrices with one row per cell.

- varm:

  Named list of numeric matrices with one row per gene.

- obsp:

  Named list of CSR matrices, each a list with `indptr`, `indices`
  (0-indexed) and `data`.

- uns_json:

  String. JSON object written to `uns`.

- chunk_size:

  Integer. Number of cells per streaming batch.

## Value

Invisible `NULL`.
