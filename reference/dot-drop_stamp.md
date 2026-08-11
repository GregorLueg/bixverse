# Strip the provenance stamp from an artefact

Used by the getters on the way out. The stamp is an internal consistency
device and has no business leaking into user code, where it would show
up in [`str()`](https://rdrr.io/r/utils/str.html) and make
[`identical()`](https://rdrr.io/r/base/identical.html) on two
numerically equal matrices `FALSE`.

## Usage

``` r
.drop_stamp(obj)
```

## Arguments

- obj:

  The artefact payload.

## Value

The payload without the stamp attribute.
