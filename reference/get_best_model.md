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

## Examples

``` r
# pull the selected fit out of a sweep without refitting
set.seed(42L)
corpus <- matrix(rbinom(200L * 40L, 1L, 0.05), nrow = 200L, ncol = 40L)
corpus[1:100, 1:10] <- rbinom(1000L, 1L, 0.6)
corpus[101:200, 11:20] <- rbinom(1000L, 1L, 0.6)
colnames(corpus) <- sprintf("term_%02d", 1:40)
sweep_res <- lda_k_sweep(corpus > 0, k_range = 5:7, .verbose = FALSE)
dim(get_best_model(sweep_res))
#> [1] 200  40   6
```
