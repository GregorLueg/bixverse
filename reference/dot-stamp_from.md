# Pull the parent artefact names out of a setter's dots

The setter generics all carry `...`, so producers opt into provenance by
passing `from = `. Producers that do not get a root stamp and keep
working, which is what makes this backwards compatible with downstream
packages calling the setters positionally.

## Usage

``` r
.stamp_from(...)
```

## Arguments

- ...:

  The setter's dots.

## Value

Character vector of parent artefact names, possibly empty.
