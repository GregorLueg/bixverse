# Helper function to generate the MAGIC imputed layer

Wraps the dense matrix
[`rs_magic_impute()`](https://gregorlueg.github.io/bixverse/reference/rs_magic_impute.md)
returns. This is the payload
[`run_magic_sc()`](https://gregorlueg.github.io/bixverse/reference/run_magic_sc.md)
writes into the cache, not a free standing result, so it is read back
with
[`get_magic()`](https://gregorlueg.github.io/bixverse/reference/get_magic.md).

## Usage

``` r
new_sc_magic(imputed, magic_params, modality)
```

## Arguments

- imputed:

  Numeric matrix of cells x genes with the imputed counts. Needs both
  dimnames set.

- magic_params:

  List. The parameters the run used, see
  [`params_sc_magic()`](https://gregorlueg.github.io/bixverse/reference/params_sc_magic.md).

- modality:

  String. The modality of the kNN graph that did the smoothing. The
  values themselves are always RNA.

## Value

Generates the `ScMagic` class.
