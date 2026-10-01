# Store the leaf coordinates of a Bonsai tree as an embedding

Writes the 2D position of every leaf, cell or metacell, into the
object's embeddings, so the usual embedding plots can use it. The tree
has to have been built over the object's current cells (or all of its
metacells), in the same order.

## Usage

``` r
set_bonsai_embedding(object, tree, name = "bonsai")
```

## Arguments

- object:

  `SingleCells` or `MetaCells` class.

- tree:

  A `BonsaiTree`, from
  [`bonsai_sc()`](https://gregorlueg.github.io/bixverse/reference/bonsai_sc.md)
  on this object.

- name:

  String. Name of the embedding.

## Value

The object with the embedding added.

## Examples

``` r
# Bonsai leaf coordinates next to the other embeddings
sc <- demo_single_cells(prepped = FALSE)
tree <- bonsai_sc(sc, .verbose = FALSE)
sc <- set_bonsai_embedding(sc, tree)
get_available_embeddings(sc)
#> [1] "bonsai"

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
