# Write an artefact payload back into a cache

Write an artefact payload back into a cache

## Usage

``` r
.sc_payload_set(cache, value, artefact, name = NULL, wnn = FALSE)
```

## Arguments

- cache:

  The cache list.

- value:

  The payload to write.

- artefact:

  String. One of `c("pca", "embedding", "knn", "snn", "magic")`.

- name:

  String. Embedding name, only used when `artefact = "embedding"`.

- wnn:

  Boolean. Whether `cache` is the `wnn` slot rather than an `ScCache`.

## Value

The cache with the payload written back.
