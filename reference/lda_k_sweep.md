# Sweep the topic count for LDA

Fits one model per requested topic count and scores each with three
model selection metrics, then combines them into a single rescaled
score. Use it to pick `k`, then pull the winning model out with
[`get_best_model()`](https://gregorlueg.github.io/bixverse/reference/get_best_model.md).

## Usage

``` r
lda_k_sweep(
  x,
  k_range,
  lda_params = params_lda(),
  top_topics_coh = 5L,
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- x:

  Documents x terms matrix. Either a `dgCMatrix`, a `dgRMatrix`, or a
  dense numeric or logical matrix. Column names are required, row names
  are generated if absent.

- k_range:

  Integer vector. The topic counts to evaluate. Every entry must be at
  least 2.

- lda_params:

  List, see
  [`params_lda()`](https://gregorlueg.github.io/bixverse/reference/params_lda.md).

- top_topics_coh:

  Integer. Number of top-scoring topics averaged into the reported
  coherence.

- seed:

  Integer. Seed for the variational initialisation.

- .verbose:

  Boolean or integer. Verbosity.

## Value

An `LdaKSweepResult`, which is a data.table with one row per `k`.

## Details

The three metrics disagree by design and are meant to be looked at
together. Arun is a symmetric KL between the singular values of the
topic-term matrix and the document length distribution (lower is
better). Cao Juan is the mean pairwise cosine similarity between topics,
so it punishes a `k` that has started splitting one signal in two (lower
is better). Mimno is UMass coherence (higher is better), averaged over
the top-scoring topics only, because coherence over all topics is
dragged down by the ones that never specialised. `combined_score` is the
min-max rescaled compromise, and `best_k` is its argmax.

One thing to know before reading `best_k`: any `k` below five is struck
out of the selection entirely, because coherence saturates on small
topic counts and would otherwise always win. That is upstream behaviour,
inherited from pycisTopic. So `best_k` can never come back below five,
however good the raw metrics of a smaller `k` look, and a sweep over
`2:6` is really a choice between five and six. The excluded rows carry
`NA` in `combined_score`, keep their metrics, and trigger a warning. If
you suspect a handful of topics, read `arun_2010` and `cao_juan_2009`
off the table yourself and pass the `k` you want to
[`get_best_model()`](https://gregorlueg.github.io/bixverse/reference/get_best_model.md).

Cost is one full fit per `k`. The corpus is built once and shared, and
each fit is already parallel over documents, so the topic counts run
sequentially rather than nesting thread pools.

## References

Arun, et al., PAKDD, 2010; Cao, et al., Neurocomputing, 2009; Mimno, et
al., EMNLP, 2011
