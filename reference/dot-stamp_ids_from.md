# Resolve parent artefact names to their current stamp ids

Producers pass parent *names* (`"pca"`, `"umap"`) rather than ids, so
they never have to reach into the stamps themselves. A name may be
modality qualified (`"adt:pca"`); otherwise `modality` is assumed. Names
that do not resolve to a stamped artefact are dropped, which is what
keeps a chain built on top of a legacy artefact from claiming a
provenance it does not have.

## Usage

``` r
.stamp_ids_from(x, from, modality = "rna")
```

## Arguments

- x:

  The object owning the cache(s).

- from:

  Character vector of parent artefact names.

- modality:

  String. Modality to resolve unqualified names against.

## Value

Character vector of stamp ids.
