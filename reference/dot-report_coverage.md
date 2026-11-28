# Warn about cells without fragments and fragments off the targets

Warn about cells without fragments and fragments off the targets

## Usage

``` r
.report_coverage(res, n_cells, what)
```

## Arguments

- res:

  List. A Rust result with `cells_missing`, `off_target_fraction` and
  `off_target_chroms`.

- n_cells:

  Integer. Number of requested cells.

- what:

  String. What the targets are, for the message.

## Value

Invisible `NULL`.
