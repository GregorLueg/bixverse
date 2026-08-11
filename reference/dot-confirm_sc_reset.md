# Ask before wiping a populated cache

Returns `TRUE` when the reset should go ahead. Skips the prompt entirely
when there is nothing to destroy, and refuses to prompt in a
non-interactive session, where
[`readline()`](https://rdrr.io/r/base/readline.html) returns `""`
immediately and would silently pick a branch nobody chose.

## Usage

``` r
.confirm_sc_reset(x, force)
```

## Arguments

- x:

  The object about to be reset.

- force:

  Boolean. Whether the caller already opted in.

## Value

Boolean. Whether to proceed.
