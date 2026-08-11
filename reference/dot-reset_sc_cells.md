# Restore every cell and clear the gene selection

Sets `cells_to_keep_idx` back to `NULL` rather than materialising the
full index vector:
[`get_cells_to_keep()`](https://gregorlueg.github.io/bixverse/reference/get_cells_to_keep.md)
already falls back to all cells on an empty slot, so this is the genuine
pristine state and costs nothing at a million cells. The HVG selection
goes too, since
[`find_hvg_sc()`](https://gregorlueg.github.io/bixverse/reference/find_hvg_sc.md)
reads the cell set and a selection made on a subset does not describe
the whole.

## Usage

``` r
.reset_sc_cells(x)
```

## Arguments

- x:

  `SingleCells` or `SingleCellsMultiModal` class.

## Value

The object with its map and DuckDB reset.
