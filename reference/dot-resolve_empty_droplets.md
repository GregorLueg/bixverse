# Resolve the empty droplet mask

Either reads the mask out of obs or hands the library sizes to Rust to
infer it. Kept separate from
[`cellsweep_sc()`](https://gregorlueg.github.io/bixverse/reference/cellsweep_sc.md)
so the resolution is testable on an obs table alone.

## Usage

``` r
.resolve_empty_droplets(obs, empty_params, .verbose = TRUE)
```

## Arguments

- obs:

  data.table. The unfiltered obs table, with `lib_size`.

- empty_params:

  List. See
  [`params_sc_empty_droplets()`](https://gregorlueg.github.io/bixverse/reference/params_sc_empty_droplets.md).
  Required: there is no safe default, since the recommended
  `method = "supplied"` needs the name of the obs column holding the
  mask.

- .verbose:

  Logical. Controls verbosity.

## Value

Logical vector, `TRUE` where the barcode is an empty droplet.
