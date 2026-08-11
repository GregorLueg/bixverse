# Stale artefacts rendered for a `print` method

Stale artefacts rendered for a `print` method

## Usage

``` r
.print_stale_str(x)
```

## Arguments

- x:

  The object owning the cache(s).

## Value

String. Comma separated labels, or `"none"`. Never signals: a `print`
method that can fail is worse than one that says "unknown".
