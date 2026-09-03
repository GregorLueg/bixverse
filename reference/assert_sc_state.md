# Assert that cached artefacts still match the object's state

The hard tier. Called at the entry of anything that hands cached indices
to Rust or writes a derived result back into the object, where a stale
artefact means either a crash in the binding layer or silently
mis-aligned biology.

## Usage

``` r
assert_sc_state(
  x,
  artefacts,
  modality = "rna",
  .var.name = checkmate::vname(x),
  add = NULL
)
```

## Arguments

- x:

  `SingleCells`, `SingleCellsSubset`, `MetaCells` or
  `SingleCellsMultiModal` class.

- artefacts:

  Character vector. Artefacts to check, either an artefact kind
  (`"pca"`, `"knn"`, `"snn"`, `"magic"`) or an embedding name
  (`"umap"`). May be modality qualified (`"adt:pca"`).

- modality:

  String. Modality to resolve unqualified names against. One of
  `c("rna", "adt", "atac", "wnn")`.

- .var.name:

  Name of the checked object to print in assertions. Defaults to the
  heuristic implemented in checkmate.

- add:

  Collection to store assertion messages. See
  [`checkmate::makeAssertCollection()`](https://mllg.github.io/checkmate/reference/AssertCollection.html).

## Value

Invisibly `x` when the assertion holds, otherwise an error.
