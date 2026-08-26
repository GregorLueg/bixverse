# Get the selected model from an LDA topic count sweep

Returns the fit at `best_k`, or at a `k` you name, without refitting.
The sweep keeps every model it fitted.

## Usage

``` r
get_best_model(x, k = NULL)

# S3 method for class 'LdaKSweepResult'
get_best_model(x, k = NULL)
```

## Arguments

- x:

  `LdaKSweepResult` object.

- k:

  Optional integer. The topic count to extract. If `NULL`, uses the
  `best_k` the sweep selected.

## Value

An `LdaResult`.

## Details

`best_k` is never below five, see
[`lda_k_sweep()`](https://gregorlueg.github.io/bixverse/reference/lda_k_sweep.md).
Pass `k` explicitly if the raw metrics point somewhere the selection
could not go.
