# Run Bonsai with a given Rust entry point

Everything
[`bonsai_sc()`](https://gregorlueg.github.io/bixverse/reference/bonsai_sc.md)
does around the Rust call: candidate genes, cells, the `BonsaiTree` and
the total timing. The Rust entry point is an argument so `bixverse.gpu`
can hand in its GPU Sanity one and get back the identical class.

## Usage

``` r
.bonsai_sc_run(object, hvg, bonsai_params, runner, .verbose)
```

## Arguments

- object:

  `SingleCells` class.

- hvg:

  Optional integer. 1-indexed candidate genes, `NULL` for all.

- bonsai_params:

  List. See
  [`params_sc_bonsai()`](https://gregorlueg.github.io/bixverse/reference/params_sc_bonsai.md).

- runner:

  Function with the signature of
  [`rs_sc_bonsai()`](https://gregorlueg.github.io/bixverse/reference/rs_sc_bonsai.md).

- .verbose:

  Boolean or integer. Controls verbosity.

## Value

A `BonsaiTree`, see
[`bonsai_sc()`](https://gregorlueg.github.io/bixverse/reference/bonsai_sc.md).
