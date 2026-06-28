# Identify HVGs without mutating object state

Like
[`find_hvg_sc()`](https://gregorlueg.github.io/bixverse/reference/find_hvg_sc.md)
but does not mutate `object`. Returns a data.table with per-gene HVG
statistics plus `is_hvg`/`hvg_rank` for the top `hvg_no` genes. Useful
for computing HVGs on a subset of cells (e.g. a specific cell type) for
downstream methods like NMF, without overwriting the HVGs stored on the
object.

## Usage

``` r
get_hvg_data_sc(
  object,
  cell_ids = NULL,
  hvg_no = 3000L,
  hvg_params = params_sc_hvg(),
  streaming = NULL,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` or `MetaCells` class.

- cell_ids:

  Optional character. Cell ids (or meta cell ids) to restrict the HVG
  calculation to. If `NULL`, uses
  [`get_cells_to_keep()`](https://gregorlueg.github.io/bixverse/reference/get_cells_to_keep.md)
  for `SingleCells` and all meta cells for `MetaCells`.

- hvg_no:

  Integer. Number of top HVGs to flag. Defaults to `3000L`.

- hvg_params:

  List, see
  [`params_sc_hvg()`](https://gregorlueg.github.io/bixverse/reference/params_sc_hvg.md).

- streaming:

  Optional Boolean. Stream the data. Ignored for `MetaCells`.

- .verbose:

  Boolean or integer. Verbosity.

## Value

data.table with `gene_idx`, `gene_id`, the HVG statistics returned by
the Rust HVG function, an `is_hvg` boolean and an `hvg_rank` integer
(`NA` for non-HVGs).
