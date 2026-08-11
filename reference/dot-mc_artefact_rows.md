# Resolve meta cell memberships onto the row space of a source artefact

`original_cell_idx` holds 1-indexed positions in the *full, unfiltered*
obs table of the source object. Embeddings, kNN graphs and diffusion
maps have one row per QC-passing cell, in `cells_to_keep` order (the
invariant the cache stamps rely on). Indexing such an artefact with the
memberships directly is silently off by the filtering the moment QC
dropped a cell, which is why the source's `cells_to_keep` is recorded on
the object at generation time. This resolves the one space into the
other.

## Usage

``` r
.mc_artefact_rows(object, n_rows)
```

## Arguments

- object:

  `MetaCells` class.

- n_rows:

  Integer. Number of rows the artefact has, checked against the recorded
  cell mapping.

## Value

A list of integer vectors, one per meta cell, holding 1-indexed row
positions in the artefact.
