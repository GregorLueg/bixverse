# lda --------------------------------------------------------------------------

## helpers ---------------------------------------------------------------------

#' Coerce a document-term matrix into the sparse list Rust wants
#'
#' @description
#' Accepts a `dgCMatrix`, a `dgRMatrix` or a dense numeric/logical matrix.
#' Logical input is what [bixverse::binarise_regulon_activity()] hands back, so
#' it is promoted rather than rejected.
#'
#' @param x The document-term matrix.
#'
#' @returns A list with `data`, `indices`, `indptr`, `cs_type`, `nrow` and
#' `ncol`, plus the row and column names as attributes `doc_ids` and
#' `term_ids`.
#'
#' @keywords internal
.lda_matrix_to_list <- function(x) {
  is_sparse <- checkmate::testClass(x, "dgCMatrix") ||
    checkmate::testClass(x, "dgRMatrix")

  if (!is_sparse) {
    # logical is what binarise_regulon_activity() hands back
    checkmate::assertMatrix(
      x,
      mode = if (is.logical(x)) "logical" else "numeric",
      min.rows = 2,
      min.cols = 2
    )
    storage.mode(x) <- "double"
    x <- methods::as(methods::as(x, "generalMatrix"), "CsparseMatrix")
  }

  doc_ids <- rownames(x)
  term_ids <- colnames(x)

  if (is.null(term_ids)) {
    stop(paste(
      "The document-term matrix needs column names. They are the terms, and a",
      "topic without term names is not interpretable."
    ))
  }
  if (is.null(doc_ids)) {
    doc_ids <- sprintf("doc_%i", seq_len(nrow(x)))
  }

  first_bad <- which(!is.finite(x@x) | x@x < 0)[1]
  if (!is.na(first_bad)) {
    # CSR stores the column per non-zero outright; CSC has to walk the column
    # pointer, which is what `@p` is there and is *not* what it is in CSR
    term <- if (inherits(x, "dgRMatrix")) {
      x@j[first_bad] + 1L
    } else {
      findInterval(first_bad, x@p, left.open = TRUE)
    }
    stop(sprintf(
      paste(
        "The document-term matrix holds negative or non-finite counts (first",
        "offending term: `%s`). LDA needs non-negative counts."
      ),
      term_ids[term]
    ))
  }

  res <- sparse_mat_to_list(x)
  attr(res, "doc_ids") <- doc_ids
  attr(res, "term_ids") <- term_ids

  res
}

#' Validate a k range for an LDA sweep
#'
#' @param k_range Integer vector. The topic counts to evaluate.
#'
#' @returns The validated integer vector.
#'
#' @keywords internal
.assert_lda_k_range <- function(k_range) {
  checkmate::assertIntegerish(
    k_range,
    lower = 2L,
    min.len = 1L,
    any.missing = FALSE,
    unique = TRUE
  )
  as.integer(k_range)
}

## main ------------------------------------------------------------------------

#' Fit a latent Dirichlet allocation model
#'
#' @description
#' Topic model over a documents x terms count matrix, fitted by variational
#' Bayes. The intended input is a binary matrix: this is the model behind
#' cisTopic, where binarised cells x regions scATAC becomes cells x topics for
#' clustering and topics x regions for region set discovery. Nothing here is
#' ATAC-specific, so binarised regulon activity from
#' [bixverse::binarise_regulon_activity()] works the same way, as does any
#' count matrix.
#'
#' @details
#' Binarising first is what makes a topic model the right tool. LDA on raw
#' single cell counts is dominated by library size, which is precisely the
#' effect cisTopic sidesteps.
#'
#' Use [bixverse::lda_k_sweep()] if you do not already know `k`.
#'
#' @param x Documents x terms matrix. Either a `dgCMatrix`, a `dgRMatrix`, or
#' a dense numeric or logical matrix. Column names are required, row names are
#' generated if absent.
#' @param k Integer. Number of topics.
#' @param lda_params List, see [bixverse::params_lda()].
#' @param seed Integer. Seed for the variational initialisation.
#' @param .verbose Boolean or integer. Verbosity.
#'
#' @returns An `LdaResult` object.
#'
#' @references Hoffman, Blei and Bach, NIPS, 2010; Bravo Gonzalez-Blas, et al.,
#' Nat Methods, 2019
#'
#' @export
run_lda <- function(
  x,
  k,
  lda_params = params_lda(),
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::qassert(k, "I1[2,)")
  assertLdaParams(lda_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  sparse_data <- .lda_matrix_to_list(x)

  lda_res <- rs_lda(
    sparse_data = sparse_data,
    k = k,
    lda_params = lda_params,
    seed = seed,
    verbose = parse_verbosity(.verbose)
  )

  new_lda_result(
    lda_res = lda_res,
    doc_ids = attr(sparse_data, "doc_ids"),
    term_ids = attr(sparse_data, "term_ids"),
    params = c(lda_params, list(k = k, seed = seed))
  )
}

#' Sweep the topic count for LDA
#'
#' @description
#' Fits one model per requested topic count and scores each with three model
#' selection metrics, then combines them into a single rescaled score. Use it
#' to pick `k`, then pull the winning model out with
#' [bixverse::get_best_model()].
#'
#' @details
#' The three metrics disagree by design and are meant to be looked at together.
#' Arun is a symmetric KL between the singular values of the topic-term matrix
#' and the document length distribution (lower is better). Cao Juan is the mean
#' pairwise cosine similarity between topics, so it punishes a `k` that has
#' started splitting one signal in two (lower is better). Mimno is UMass
#' coherence (higher is better), averaged over the top-scoring topics only,
#' because coherence over all topics is dragged down by the ones that never
#' specialised. `combined_score` is the min-max rescaled compromise, and
#' `best_k` is its argmax.
#'
#' One thing to know before reading `best_k`: any `k` below five is struck out
#' of the selection entirely, because coherence saturates on small topic counts
#' and would otherwise always win. That is upstream behaviour, inherited from
#' pycisTopic. So `best_k` can never come back below five, however good the raw
#' metrics of a smaller `k` look, and a sweep over `2:6` is really a choice
#' between five and six. The excluded rows carry `NA` in `combined_score`, keep
#' their metrics, and trigger a warning. If you suspect a handful of topics,
#' read `arun_2010` and `cao_juan_2009` off the table yourself and pass the `k`
#' you want to [bixverse::get_best_model()].
#'
#' Cost is one full fit per `k`. The corpus is built once and shared, and each
#' fit is already parallel over documents, so the topic counts run sequentially
#' rather than nesting thread pools.
#'
#' @inheritParams run_lda
#'
#' @param k_range Integer vector. The topic counts to evaluate. Every entry
#' must be at least 2.
#' @param top_topics_coh Integer. Number of top-scoring topics averaged into
#' the reported coherence.
#'
#' @returns An `LdaKSweepResult`, which is a data.table with one row per `k`.
#'
#' @references Arun, et al., PAKDD, 2010; Cao, et al., Neurocomputing, 2009;
#' Mimno, et al., EMNLP, 2011
#'
#' @export
lda_k_sweep <- function(
  x,
  k_range,
  lda_params = params_lda(),
  top_topics_coh = 5L,
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  k_range <- .assert_lda_k_range(k_range)
  assertLdaParams(lda_params)
  checkmate::qassert(top_topics_coh, "I1[1,)")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  sparse_data <- .lda_matrix_to_list(x)

  sweep_res <- rs_lda_k_sweep(
    sparse_data = sparse_data,
    k_range = k_range,
    lda_params = lda_params,
    top_topics_coh = as.integer(top_topics_coh),
    seed = seed,
    verbose = parse_verbosity(.verbose)
  )

  new_lda_k_sweep_result(
    sweep_res = sweep_res,
    doc_ids = attr(sparse_data, "doc_ids"),
    term_ids = attr(sparse_data, "term_ids"),
    params = c(lda_params, list(seed = seed, top_topics_coh = top_topics_coh))
  )
}
