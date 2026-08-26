# Fit a latent Dirichlet allocation model to a document-term matrix

**\[experimental\]** Variational Bayes, following Hoffman, et al. The
matrix is documents x terms; for a cisTopic-style run that is a
binarised cells x regions (or cells x regulons) matrix, but any count
matrix works.

## Usage

``` r
rs_lda(sparse_data, k, lda_params, seed, verbose)
```

## Arguments

- sparse_data:

  A named list that needs to have `data`, `indptr`, `indices`, `nrow`,
  `ncol` and `cs_type`. Documents x terms.

- k:

  Integer. Number of topics.

- lda_params:

  List. The LDA parameters, see
  [`params_lda()`](https://gregorlueg.github.io/bixverse/reference/params_lda.md).

- seed:

  Integer. Seed for the variational initialisation.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following items

- cell_topic - Topic proportions per document, `k x n_documents`.

- topic_region - Term probabilities per topic, `n_terms x k`.

- bound - Final variational bound (ELBO). Higher is better.

- perplexity - Per-token perplexity. Lower is better.

- n_iter - Outer iterations run.

- converged - Whether the relative bound change fell below `tol`.

## References

Hoffman, Blei and Bach, NIPS, 2010
