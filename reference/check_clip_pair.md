# Check that a clipping range is whole

Rust falls back to its own default when only one end is given, so half a
range silently becomes no range at all.

## Usage

``` r
check_clip_pair(x, label)
```

## Arguments

- x:

  The parameter list.

- label:

  Short label used in the error message.

## Value

`TRUE` if the check was successful, otherwise an error message.
