# Attach a provenance stamp to an artefact

Attach a provenance stamp to an artefact

## Usage

``` r
.set_stamp(obj, stamp)
```

## Arguments

- obj:

  The artefact payload (matrix, kNN object, igraph).

- stamp:

  List. The stamp as built by
  [`.new_sc_stamp()`](https://gregorlueg.github.io/bixverse/reference/dot-new_sc_stamp.md).

## Value

The payload with the stamp attached as an attribute.
