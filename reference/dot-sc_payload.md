# Read an artefact payload straight out of a cache

Raw slot access with no freshness checking and no stamp stripping. This
is what the provenance machinery reads; user facing code goes through
the getters.

The `wnn` pseudo modality stores its embeddings under `embeddings`,
where a real `ScCache` uses `other_embeddings`, hence the flag.

## Usage

``` r
.sc_payload(cache, artefact, name = NULL, wnn = FALSE)
```

## Arguments

- cache:

  The cache list.

- artefact:

  String. One of
  `c("pca", "embedding", "knn", "connectivities", "magic", "residuals")`.

- name:

  String. Embedding name, only used when `artefact = "embedding"`.

- wnn:

  Boolean. Whether `cache` is the `wnn` slot rather than an `ScCache`.

## Value

The payload, or `NULL`.
