# Run Bonsai over metacells with a given Rust entry point

The `MetaCells` counterpart of
[`.bonsai_sc_run()`](https://gregorlueg.github.io/bixverse/reference/dot-bonsai_sc_run.md).
The raw counts go to Rust from memory, and every metacell is a leaf.

## Usage

``` r
.bonsai_mc_run(object, hvg, bonsai_params, runner, .verbose)
```

## Arguments

- object:

  `MetaCells` class.

- hvg:

  Optional integer. 1-indexed candidate genes, `NULL` for all.

- bonsai_params:

  List. See
  [`params_sc_bonsai()`](https://gregorlueg.github.io/bixverse/reference/params_sc_bonsai.md).

- runner:

  Function with the signature of
  [`rs_mc_bonsai()`](https://gregorlueg.github.io/bixverse/reference/rs_mc_bonsai.md).

- .verbose:

  Boolean or integer. Controls verbosity.

## Value

A `BonsaiTree`, see
[`bonsai_sc()`](https://gregorlueg.github.io/bixverse/reference/bonsai_sc.md).
