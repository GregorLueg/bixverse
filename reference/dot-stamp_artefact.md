# Stamp an artefact that has just been written into a cache

Called by the S7 setter forwarders, which are the only layer with access
to both the cache and the cell state. A no-op when the artefact is
absent.

## Usage

``` r
.stamp_artefact(x, artefact, name = NULL, modality = "rna", from = character())
```

## Arguments

- x:

  The object owning the cache.

- artefact:

  String. One of `c("pca", "embedding", "knn", "snn", "magic")`.

- name:

  String. Embedding name, only used when `artefact = "embedding"`.

- modality:

  String. The modality the artefact was written to.

- from:

  Character vector of parent artefact names.

## Value

The object with the artefact stamped.
