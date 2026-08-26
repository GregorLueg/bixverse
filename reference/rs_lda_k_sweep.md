# Fit LDA across a range of topic counts and score each fit

**\[experimental\]** Fits one model per requested topic count and
evaluates each with the Arun, Cao Juan and Mimno metrics. The corpus is
built once and shared, and each fit is already parallel over documents,
so the topic counts run sequentially rather than nesting Rayon pools.

## Usage

``` r
rs_lda_k_sweep(sparse_data, k_range, lda_params, top_topics_coh, seed, verbose)
```

## Arguments

- sparse_data:

  A named list that needs to have `data`, `indptr`, `indices`, `nrow`,
  `ncol` and `cs_type`. Documents x terms.

- k_range:

  Integer vector. Topic counts to evaluate.

- lda_params:

  List. The LDA parameters, see
  [`params_lda()`](https://gregorlueg.github.io/bixverse/reference/params_lda.md).

- top_topics_coh:

  Optional integer. Number of top-scoring topics averaged into the
  reported coherence. If `NULL`, defaults to `5L`.

- seed:

  Integer. Seed for the variational initialisation.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following items

- k - The topic counts that were tried.

- models - One fitted model per topic count, see
  [`rs_lda()`](https://gregorlueg.github.io/bixverse/reference/rs_lda.md).

- metrics - One metric list per topic count, with `arun_2010`,
  `cao_juan_2009`, `mimno_2011`, `coherence_per_topic`, `bound` and
  `perplexity`.

- combined_score - Rescaled combined score per entry. `NaN` where the
  entry was excluded from selection by the coherence topic-count floor.

- best_k - Topic count with the highest combined score.

## References

Arun, et al., PAKDD, 2010; Cao, et al., Neurocomputing, 2009; Mimno, et
al., EMNLP, 2011
