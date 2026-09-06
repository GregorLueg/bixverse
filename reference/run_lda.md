# Fit a latent Dirichlet allocation model

Topic model over a documents x terms count matrix, fitted by variational
Bayes. The intended input is a binary matrix: this is the model behind
cisTopic, where binarised cells x regions scATAC becomes cells x topics
for clustering and topics x regions for region set discovery. Nothing
here is ATAC-specific, so binarised regulon activity from
[`binarise_regulon_activity()`](https://gregorlueg.github.io/bixverse/reference/binarise_regulon_activity.md)
works the same way, as does any count matrix.

## Usage

``` r
run_lda(x, k, lda_params = params_lda(), seed = 42L, .verbose = TRUE)
```

## Arguments

- x:

  Documents x terms matrix. Either a `dgCMatrix`, a `dgRMatrix`, or a
  dense numeric or logical matrix. Column names are required, row names
  are generated if absent.

- k:

  Integer. Number of topics.

- lda_params:

  List, see
  [`params_lda()`](https://gregorlueg.github.io/bixverse/reference/params_lda.md).

- seed:

  Integer. Seed for the variational initialisation.

- .verbose:

  Boolean or integer. Verbosity.

## Value

An `LdaResult` object.

## Details

Binarising first is what makes a topic model the right tool. LDA on raw
single cell counts is dominated by library size, which is precisely the
effect cisTopic sidesteps.

Use
[`lda_k_sweep()`](https://gregorlueg.github.io/bixverse/reference/lda_k_sweep.md)
if you do not already know `k`.

## References

Hoffman, Blei and Bach, NIPS, 2010; Bravo Gonzalez-Blas, et al., Nat
Methods, 2019

## Examples

``` r
# two planted term blocks recovered as two topics
set.seed(42L)
corpus <- matrix(rbinom(200L * 40L, 1L, 0.05), nrow = 200L, ncol = 40L)
corpus[1:100, 1:10] <- rbinom(1000L, 1L, 0.6)
corpus[101:200, 11:20] <- rbinom(1000L, 1L, 0.6)
colnames(corpus) <- sprintf("term_%02d", 1:40)
lda_res <- run_lda(corpus > 0, k = 2L, .verbose = FALSE)
lda_res
#> LdaResult (latent Dirichlet allocation)
#>   Documents:        200
#>   Terms:            40
#>   Topics:           2
#>   Bound (ELBO):     -5137.97
#>   Perplexity:       31.6666
#>   Iterations:       60
#> 
```
