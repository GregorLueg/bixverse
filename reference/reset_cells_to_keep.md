# Reset the cells to keep

Restores every cell found in the binary count file and wipes the cache,
taking the object back to a pristine state. Filtering only ever flips a
`to_keep` flag in the DuckDB, so nothing was deleted and nothing is lost
by resetting.

Wiping the cache is not optional: a PCA or a kNN computed on a filtered
subset does not describe the full cell set, and keeping it would
recreate exactly the mismatch this guards against. With `force = FALSE`
you are asked to confirm before that happens.

## Usage

``` r
reset_cells_to_keep(object, force = FALSE)

## S7 method for class <bixverse::MetaCells>
reset_cells_to_keep(object, force = FALSE)

## S7 method for class <bixverse::SingleCells>
reset_cells_to_keep(object, force = FALSE)

## S7 method for class <bixverse::SingleCellsMultiModal>
reset_cells_to_keep(object, force = FALSE)

## S7 method for class <bixverse::SingleCellsSubset>
reset_cells_to_keep(object, force = FALSE)
```

## Arguments

- object:

  `SingleCells` or `SingleCellsMultiModal` class.

- force:

  Boolean. Skip the confirmation prompt. Defaults to `FALSE`, in which
  case an interactive session asks before wiping the cache and a
  non-interactive one errors, because there is no one there to ask.

## Value

The object with every cell restored and an empty cache. Unchanged if the
confirmation was declined.
