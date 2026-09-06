# Save memory-bound data to disk

Helper function that stores the memory-bound data to disk for
checkpointing or when you close the session for quick recovery of prior
work. You have the option to save as `".rds"` or `".qs2"` (you need to
have the package `"qs2"` installed for this option!).

## Usage

``` r
save_sc_exp_to_disk(object, type = c("qs2", "rds"))
```

## Arguments

- object:

  `SingleCells`, `MetaCells` or `SingleCellsMultiModal` class.

- type:

  String. One of `c("qs2", "rds")`. Defines which binary format to use.
  Will default to `"qs2"` for speed.

## Value

`NULL`, invisibly. Called for the side effect of writing the in-memory
maps and caches next to the counts. It does not return the object, so do
not assign the result.

## Examples

``` r
# checkpoint the in-memory map and cache next to the counts
sc <- demo_single_cells()
save_sc_exp_to_disk(sc, type = "rds")
list.files(sc@dir_data)
#> [1] "counts_cells.bin" "counts_genes.bin" "memory.rds"       "sc_duckdb.db"    

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
