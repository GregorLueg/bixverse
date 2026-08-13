# Resolve the shared inputs of the meta cell HotSpot methods

Both HotSpot entry points need the same four things resolved off the
object before anything can go to Rust. Kept in one place so the two
methods cannot drift apart.

## Usage

``` r
.prep_mc_hotspot(
  object,
  embd_to_use,
  use_knn,
  no_embd_to_use,
  cells_to_take,
  genes_to_take
)
```

## Arguments

- object:

  `MetaCells` class.

- embd_to_use:

  String. The embedding to use.

- use_knn:

  Boolean. Shall the cached kNN graph be used. Forced to `FALSE` when
  `cells_to_take` is provided, since the cached graph covers all meta
  cells.

- no_embd_to_use:

  Optional integer. Number of embedding dimensions to use. If `NULL` all
  will be used.

- cells_to_take:

  Optional string vector. Meta cell identifiers. If `NULL` all meta
  cells are used.

- genes_to_take:

  Optional string vector. If `NULL` all genes are used.

## Value

A list with the following elements:

- embd - The embedding matrix, subset to `cells_to_take`.

- cells_to_take - The resolved meta cell identifiers.

- genes_to_take - The resolved gene identifiers.

- knn_data - The cached kNN object or `NULL`.
