# Full state report for an object's caches

Full state report for an object's caches

## Usage

``` r
.sc_state_report(x)
```

## Arguments

- x:

  The object owning the cache(s).

## Value

A `data.table` with one row per present artefact and the columns

- modality - The modality the artefact lives in.

- artefact - One of `pca`, `embedding`, `knn`, `snn`, `magic`.

- name - The embedding name, `NA` for the others.

- stamped - Whether the artefact carries a provenance stamp.

- stale - Whether it disagrees with the current state.

- reason - Why it is stale, `NA` otherwise.

- id - The artefact's stamp id.

- from - List column of the parent ids.
