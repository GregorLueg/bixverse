# Why an artefact is stale, if it is

Walks the provenance chain. An artefact is stale when the cell set moved
under it, when a parent it was built on no longer exists (i.e. was
recomputed and minted a new id), or when a parent is itself stale. The
last case is why this recurses: without it, re-running the PCA marks the
kNN stale but the sNN still resolves its untouched parent and passes.

## Usage

``` r
.sc_stale_reason(label, artefacts, cells_hash, seen = character())
```

## Arguments

- label:

  String. The artefact label to evaluate.

- artefacts:

  List as returned by
  [`.sc_artefacts()`](https://gregorlueg.github.io/bixverse/reference/dot-sc_artefacts.md).

- cells_hash:

  String. The object's current cell state hash.

- seen:

  Character vector. Labels already on the recursion stack.

## Value

`NA_character_` when fresh (or unstamped), otherwise the reason.
