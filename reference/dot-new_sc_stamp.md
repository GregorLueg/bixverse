# Build a provenance stamp

Build a provenance stamp

## Usage

``` r
.new_sc_stamp(cells, from = character())
```

## Arguments

- cells:

  String. Hash of the cell state at write time.

- from:

  Character vector. Ids of the artefacts this one derives from. Empty
  for a root artefact such as the PCA.

## Value

A list with `cells`, `id` and `from`.
