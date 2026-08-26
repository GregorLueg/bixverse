# lda --------------------------------------------------------------------------

set.seed(123L)

## planted corpus --------------------------------------------------------------

# two disjoint term blocks firing in two disjoint document blocks, on a thin
# background. Any working topic model has to put the two blocks in two topics.
n_docs <- 400L
n_terms <- 60L
block <- 20L

corpus <- matrix(
  rbinom(n_docs * n_terms, 1L, 0.02),
  nrow = n_docs,
  ncol = n_terms
)
corpus[1:200, 1:block] <- rbinom(200L * block, 1L, 0.6)
corpus[201:400, (block + 1L):(2L * block)] <- rbinom(200L * block, 1L, 0.6)
colnames(corpus) <- sprintf("term_%02d", seq_len(n_terms))
rownames(corpus) <- sprintf("doc_%03d", seq_len(n_docs))

block_a <- colnames(corpus)[1:block]
block_b <- colnames(corpus)[(block + 1L):(2L * block)]

lda_res <- run_lda(corpus > 0, k = 2L, .verbose = FALSE)

expect_inherits(
  lda_res,
  "LdaResult",
  info = "run_lda() returns an LdaResult"
)

## planted topic recovery ------------------------------------------------------

top_terms <- get_top_terms(lda_res, n = block)

# ground truth decides which topic is which; everything below is a statement
# about how cleanly the two blocks separated
hits_a <- vapply(
  c("topic_01", "topic_02"),
  \(tp) sum(top_terms$term[top_terms$topic == tp] %in% block_a),
  integer(1)
)

expect_equal(
  current = as.integer(sort(hits_a, decreasing = TRUE)),
  target = c(block, 0L),
  info = paste(
    "one topic's top terms are entirely block A and the other's hold none",
    "of it, i.e. the planted blocks landed in separate topics"
  )
)

pure_topic <- names(which.max(hits_a))
other_topic <- setdiff(c("topic_01", "topic_02"), pure_topic)

expect_equal(
  current = sum(top_terms$term[top_terms$topic == other_topic] %in% block_b),
  target = block,
  info = "the other topic's top terms are entirely block B"
)

# documents 1:200 carry block A, so they should load on the block A topic
doc_topic <- as.matrix(lda_res, "doc_topic")

expect_true(
  current = all(doc_topic[1:200, pure_topic] > 0.5),
  info = "every block A document puts most of its mass on the block A topic"
)

expect_true(
  current = all(doc_topic[201:400, other_topic] > 0.5),
  info = "every block B document puts most of its mass on the block B topic"
)

## shapes, names and simplex ---------------------------------------------------

expect_equal(
  current = dim(lda_res),
  target = c(n_docs, n_terms, 2L),
  info = "dim() gives documents, terms and topics"
)

expect_equal(
  current = dimnames(doc_topic),
  target = list(rownames(corpus), c("topic_01", "topic_02")),
  info = "doc_topic carries the document ids and topic names"
)

expect_equal(
  current = dimnames(as.matrix(lda_res, "term_topic")),
  target = list(colnames(corpus), c("topic_01", "topic_02")),
  info = "term_topic carries the term ids and topic names"
)

expect_true(
  current = all(abs(rowSums(doc_topic) - 1) < 1e-6),
  info = "each document's topic proportions sum to one"
)

expect_true(
  current = all(abs(colSums(as.matrix(lda_res, "term_topic")) - 1) < 1e-6),
  info = "each topic's term probabilities sum to one"
)

## reproducibility -------------------------------------------------------------

expect_identical(
  current = run_lda(corpus > 0, k = 2L, .verbose = FALSE)$doc_topic,
  target = lda_res$doc_topic,
  info = "the same seed gives an identical fit"
)

expect_false(
  current = identical(
    run_lda(corpus > 0, k = 2L, seed = 99L, .verbose = FALSE)$doc_topic,
    lda_res$doc_topic
  ),
  info = "a different seed gives a different initialisation"
)

## input coercion --------------------------------------------------------------

corpus_num <- corpus * 1
corpus_csc <- methods::as(
  methods::as(corpus_num, "generalMatrix"),
  "CsparseMatrix"
)

expect_equal(
  current = run_lda(corpus_csc, k = 2L, .verbose = FALSE)$doc_topic,
  target = lda_res$doc_topic,
  info = "a dgCMatrix gives the same fit as the dense logical matrix"
)

expect_equal(
  current = run_lda(
    methods::as(corpus_csc, "RsparseMatrix"),
    k = 2L,
    .verbose = FALSE
  )$doc_topic,
  target = lda_res$doc_topic,
  info = "a dgRMatrix gives the same fit, i.e. both orientations are handled"
)

## input validation ------------------------------------------------------------

expect_error(
  run_lda(corpus_num, k = 1L, .verbose = FALSE),
  pattern = "k",
  info = "a single topic is rejected"
)

negative <- corpus_num
negative[1L, 3L] <- -1
expect_error(
  run_lda(negative, k = 2L, .verbose = FALSE),
  pattern = "term_03",
  info = "negative counts are caught in R and name the offending term"
)

# CSR stores the column index per non-zero, CSC has to walk the column pointer.
# Getting that wrong names a plausible but incorrect term, so check both.
negative_csr <- methods::as(
  methods::as(methods::as(negative, "generalMatrix"), "CsparseMatrix"),
  "RsparseMatrix"
)
expect_error(
  run_lda(negative_csr, k = 2L, .verbose = FALSE),
  pattern = "term_03",
  info = "the offending term is named correctly for CSR input too"
)

unnamed <- corpus_num
colnames(unnamed) <- NULL
expect_error(
  run_lda(unnamed, k = 2L, .verbose = FALSE),
  pattern = "needs column names",
  info = "a matrix without term names is rejected"
)

no_rownames <- corpus_num
rownames(no_rownames) <- NULL
expect_equal(
  current = rownames(run_lda(no_rownames, k = 2L, .verbose = FALSE)$doc_topic),
  target = sprintf("doc_%i", seq_len(n_docs)),
  info = "missing document names are generated rather than rejected"
)

expect_error(
  run_lda(corpus_num, k = 2L, lda_params = list(alpha = 1)),
  pattern = "Names must include",
  info = "an incomplete parameter list is rejected"
)

## parameter wrapper -----------------------------------------------------------

expect_error(
  params_lda(alpha = -1),
  pattern = "alpha",
  info = "params_lda() rejects a non-positive alpha"
)

expect_error(
  params_lda(learning = "gibbs"),
  info = "params_lda() rejects an unknown learning strategy"
)

expect_true(
  current = isTRUE(bixverse:::checkLdaParams(params_lda())),
  info = "the defaults pass their own assertion"
)

broken <- params_lda()
broken$check_every <- -1L
expect_true(
  current = is.character(bixverse:::checkLdaParams(broken)),
  info = "checkLdaParams() reports a message for a broken element"
)

broken_choice <- params_lda()
broken_choice$learning <- "gibbs"
expect_true(
  current = is.character(bixverse:::checkLdaParams(broken_choice)),
  info = "checkLdaParams() reports a message for an unknown learning strategy"
)

## path dispatch ---------------------------------------------------------------

online_res <- run_lda(
  corpus > 0,
  k = 2L,
  lda_params = params_lda(
    learning = "online",
    n_epochs = 5L,
    batch_size = 64L
  ),
  .verbose = FALSE
)

expect_inherits(
  online_res,
  "LdaResult",
  info = "the online learning path runs"
)

expect_true(
  current = is.finite(online_res$perplexity) && online_res$perplexity > 0,
  info = "the online fit reports a usable perplexity"
)

## k sweep ---------------------------------------------------------------------

k_range <- 5:8

# the coherence topic-count floor is 5, so this range is entirely eligible and
# the sweep runs without the exclusion warning
sweep_res <- lda_k_sweep(corpus > 0, k_range = k_range, .verbose = FALSE)

expect_inherits(
  sweep_res,
  "LdaKSweepResult",
  info = "lda_k_sweep() returns an LdaKSweepResult"
)

expect_inherits(
  sweep_res,
  "data.table",
  info = "the sweep result is also a data.table"
)

expect_equal(
  current = sweep_res$k,
  target = as.integer(k_range),
  info = "one row per requested topic count, in the order requested"
)

expect_true(
  current = all(
    c(
      "arun_2010",
      "cao_juan_2009",
      "mimno_2011",
      "bound",
      "perplexity",
      "combined_score",
      "converged"
    ) %in%
      names(sweep_res)
  ),
  info = "the sweep reports every metric"
)

expect_true(
  current = attr(sweep_res, "best_k") %in% k_range,
  info = "best_k is one of the swept topic counts"
)

best_model <- get_best_model(sweep_res)

expect_inherits(
  best_model,
  "LdaResult",
  info = "get_best_model() returns an LdaResult"
)

expect_equal(
  current = ncol(as.matrix(best_model, "doc_topic")),
  target = attr(sweep_res, "best_k"),
  info = "the returned model has best_k topics"
)

expect_equal(
  current = ncol(as.matrix(get_best_model(sweep_res, k = 6L), "doc_topic")),
  target = 6L,
  info = "an explicit k pulls that fit out instead"
)

expect_error(
  get_best_model(sweep_res, k = 42L),
  pattern = "No fit for k",
  info = "asking for a k that was never fitted is an error"
)

## coherence floor -------------------------------------------------------------

# k below 5 cannot be selected, so a sweep that straddles the floor has to say
# so rather than quietly returning a k the metrics never argued for
expect_warning(
  lda_k_sweep(corpus > 0, k_range = 3:6, .verbose = FALSE),
  pattern = "coherence topic-count floor",
  info = "a sweep straddling the floor warns about the excluded entries"
)

straddling <- suppressWarnings(
  lda_k_sweep(corpus > 0, k_range = 3:6, .verbose = FALSE)
)

expect_true(
  current = all(is.na(straddling$combined_score[straddling$k < 5L])),
  info = "entries below the floor carry NA in combined_score"
)

expect_true(
  current = all(!is.na(straddling$arun_2010)),
  info = "entries below the floor still report their raw metrics"
)

expect_true(
  current = attr(straddling, "best_k") >= 5L,
  info = "best_k never comes back below the coherence floor"
)

expect_error(
  lda_k_sweep(corpus > 0, k_range = c(1L, 3L), .verbose = FALSE),
  info = "a k range reaching below 2 is rejected"
)

## getters ---------------------------------------------------------------------

expect_equal(
  current = get_params(lda_res)$k,
  target = 2L,
  info = "get_params() dispatches on LdaResult and carries k"
)

expect_equal(
  current = get_params(sweep_res)$seed,
  target = 42L,
  info = "get_params() dispatches on LdaKSweepResult"
)

expect_equal(
  current = nrow(get_top_terms(lda_res, n = 5L)),
  target = 10L,
  info = "get_top_terms() returns n rows per topic"
)

expect_true(
  current = all(
    c("topic", "rank", "term", "probability") %in%
      names(get_top_terms(lda_res, n = 5L))
  ),
  info = "get_top_terms() has the documented columns"
)

expect_equal(
  current = nrow(get_top_terms(lda_res, n = n_terms + 10L)),
  target = 2L * n_terms,
  info = "asking for more terms than the vocabulary holds is capped, not an error"
)
