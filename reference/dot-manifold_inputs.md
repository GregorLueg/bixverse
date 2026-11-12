# Resolve the kNN and input embedding for a manifold method

Shared entry point of the `*_sc` 2D embedding methods. `"wnn"` takes the
integrated kNN graph but reads the input embedding from the RNA cache,
so the modality the embedding comes from (`cache_modality`) can differ
from the one the result is written to.

## Usage

``` r
.manifold_inputs(
  object,
  use_knn,
  embd_to_use,
  no_embd_to_use,
  modality,
  k = NULL
)
```

## Arguments

- object:

  `SingleCells`, `MetaCells` or `SingleCellsSubset` class.

- use_knn:

  Boolean. Use the kNN graph found in the object.

- embd_to_use:

  String. The embedding to feed the method. Must be available in the
  object.

- no_embd_to_use:

  Optional integer. Number of embedding dimensions to use. If `NULL` all
  will be used.

- modality:

  String. One of `c("rna", "adt", "wnn")`.

- k:

  Optional integer. Neighbour count the caller asked for. Only feeds the
  mismatch warning and the returned `k`.

## Value

A list with

- knn - The manifoldsR `NearestNeighbours`, or `NULL` if the method
  should build its own.

- embd - The input embedding, cells x dimensions.

- cache_modality - The modality the embedding was read from.

- k - The neighbour count to hand to manifoldsR: the cached graph's
  width when one is used, `k` otherwise. manifoldsR cuts a wider graph
  down to `k`, so passing the width keeps the whole cached graph.
