# Resolve the configured strictness

`bixverse.cache_check` accepts `"error"`, `"warn"` and `"none"`. Unset
means the two tier default: the getters warn,
[`assert_sc_state()`](https://gregorlueg.github.io/bixverse/reference/assert_sc_state.md)
errors. `"error"` promotes the getter warnings to errors, `"warn"` is
the explicit spelling of the default getter behaviour, and `"none"`
disables both tiers.

## Usage

``` r
.sc_check_mode()
```

## Value

String. One of `c("default", "error", "warn", "none")`.
