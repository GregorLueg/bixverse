# Check that cached artefacts still match the object's state

Returns `TRUE` when every requested artefact is fresh, otherwise a
message describing the first problem. An artefact that is absent, or
that predates provenance stamping, passes: unknown provenance is not the
same as known bad.

Named in snake case rather than the `checkXxx` camel case of the other
checkmate extensions in this package, because it validates object state
rather than a `params_*()` list and sits alongside
[`get_sc_cache_status()`](https://gregorlueg.github.io/bixverse/reference/get_sc_cache_status.md).

## Usage

``` r
check_sc_state(x, artefacts, modality = "rna")
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

## Value

`TRUE` if the check was successful, otherwise a string describing the
failure.
