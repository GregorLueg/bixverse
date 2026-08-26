# Wrapper function for the LDA parameters

Solver options for the variational Bayes latent Dirichlet allocation,
see
[`run_lda()`](https://gregorlueg.github.io/bixverse/reference/run_lda.md).

## Usage

``` r
params_lda(
  alpha = 50,
  alpha_by_topic = TRUE,
  eta = 0.1,
  eta_by_topic = FALSE,
  max_iter = 150L,
  tol = 0.001,
  inner_max_iter = 100L,
  inner_tol = 0.001,
  check_every = 10L,
  learning = c("batch", "online"),
  batch_size = 1024L,
  n_epochs = 10L
)
```

## Arguments

- alpha:

  Float. Dirichlet prior on the document-topic distributions.

- alpha_by_topic:

  Boolean. Shall `alpha` be divided by the topic count.

- eta:

  Float. Dirichlet prior on the topic-term distributions.

- eta_by_topic:

  Boolean. Shall `eta` be divided by the topic count.

- max_iter:

  Integer. Maximum outer iterations. Ignored by the online variant,
  which counts epochs instead.

- tol:

  Float. Relative change in the bound below which the solver stops.

- inner_max_iter:

  Integer. Maximum fixed-point iterations of the per-document E-step.

- inner_tol:

  Float. Relative L1 change in the variational parameters below which
  the per-document E-step stops.

- check_every:

  Integer. Iterations between bound evaluations.

- learning:

  String. One of `c("batch", "online")`.

- batch_size:

  Integer. Documents per mini-batch. Online only.

- n_epochs:

  Integer. Passes over the corpus. Online only.

## Value

A list with the LDA parameters.

## Details

The defaults follow pycisTopic, so the knobs mean the same thing on both
sides. `alpha_by_topic = TRUE` turns `alpha` into the Griffiths and
Steyvers `50 / k` heuristic that cisTopic defaults to; set it to `FALSE`
if you want `alpha` taken literally.

`learning = "batch"` sweeps every document once per iteration and is
monotone in the bound. `"online"` takes decaying steps from shuffled
mini-batches, which reaches a usable fit in far fewer passes on a large
corpus at the cost of that guarantee. `batch_size` and `n_epochs` are
only read by the online variant.

## References

Hoffman, Blei and Bach, NIPS, 2010; Bravo Gonzalez-Blas, et al., Nat
Methods, 2019
