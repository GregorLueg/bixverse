# Name a manifold embedding and write it back onto the object

Name a manifold embedding and write it back onto the object

## Usage

``` r
.store_manifold(object, embd, prefix, slot_name, modality, from)
```

## Arguments

- object:

  `SingleCells`, `MetaCells` or `SingleCellsSubset` class.

- embd:

  Numeric matrix. The embedding, cells x dimensions.

- prefix:

  String. Column name prefix, e.g. `"umap"`.

- slot_name:

  String. Name of the embedding within the object.

- modality:

  String. Modality the embedding is written to.

- from:

  Character vector. Parent artefact names for the provenance stamp, see
  [`.manifold_from()`](https://gregorlueg.github.io/bixverse/reference/dot-manifold_from.md).

## Value

The object with the embedding added.
