# Is an artefact present in a cache?

Pure presence probe with no freshness checking. For the call sites that
ask whether something is there before overwriting it, where erroring on
a stale artefact that is about to be replaced would be actively wrong.

## Usage

``` r
.sc_has_artefact(x, artefact, name = NULL, modality = "rna")
```

## Arguments

- x:

  The object owning the cache.

- artefact:

  String. One of `c("pca", "embedding", "knn", "snn", "magic")`.

- name:

  String. Embedding name, only used when `artefact = "embedding"`.

- modality:

  String. The modality to look in.

## Value

Boolean.
