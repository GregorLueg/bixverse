# Enumerate every artefact present on an object

Enumerate every artefact present on an object

## Usage

``` r
.sc_artefacts(x)
```

## Arguments

- x:

  The object owning the cache(s).

## Value

A list of lists, each with `modality`, `artefact`, `name`, `label` and
`stamp` (which is `NULL` for an unstamped, i.e. legacy, artefact).
