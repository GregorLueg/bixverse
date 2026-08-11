# Warn (or error) when a getter hands back a stale artefact

The soft tier. Every consumer reads through the getters, so this catches
everything, but it stays a warning because several call sites use the
getters as presence probes on artefacts they are about to overwrite.

## Usage

``` r
.warn_sc_state(x, artefact, name = NULL, modality = "rna")
```

## Arguments

- x:

  The object owning the cache.

- artefact:

  String. One of `c("pca", "embedding", "knn", "snn", "magic")`.

- name:

  String. Embedding name, only used when `artefact = "embedding"`.

- modality:

  String. The modality read from.

## Value

Invisibly `TRUE`. Called for the side effect.
