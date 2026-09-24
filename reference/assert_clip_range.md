# Validate a residual clipping range

Rust falls back to its own default when only one end is supplied, so
half a range silently becomes no range at all. Catch it here instead.

## Usage

``` r
assert_clip_range(clip_min, clip_max)
```

## Arguments

- clip_min:

  Float or `NULL`. Lower bound.

- clip_max:

  Float or `NULL`. Upper bound.

## Value

Invisibly `TRUE`; called for the error.
