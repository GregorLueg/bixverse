# lda result classes -----------------------------------------------------------

## constants -------------------------------------------------------------------

# Mirrors DEFAULT_MIN_TOPICS_COH on the Rust side, which follows pycisTopic's
# min_topics_coh. Kept here only so the warnings can name the number.
LDA_MIN_TOPICS_COH <- 5L

## LdaResult -------------------------------------------------------------------

#' Generate an LdaResult instance
#'
#' @description
#' Takes the raw Rust output of [bixverse::rs_lda()] and names both matrices.
#' Rust returns the document-topic matrix as `k x n_documents` because a
#' document's topic vector is one contiguous column there; this flips it so
#' both matrices are entity-major, matching the NMF result classes.
#'
#' @param lda_res List. The raw return of [bixverse::rs_lda()].
#' @param doc_ids Character vector. The document identifiers, in row order of
#' the input matrix.
#' @param term_ids Character vector. The term identifiers, in column order of
#' the input matrix.
#' @param params List. The parameters the model was fitted with.
#'
#' @returns An object of class `LdaResult`.
#'
#' @export
#'
#' @keywords internal
new_lda_result <- function(lda_res, doc_ids, term_ids, params) {
  # checks
  checkmate::assertList(lda_res)
  checkmate::assertNames(
    names(lda_res),
    must.include = c(
      "cell_topic",
      "topic_region",
      "bound",
      "perplexity",
      "n_iter",
      "converged"
    )
  )
  checkmate::qassert(doc_ids, "S+")
  checkmate::qassert(term_ids, "S+")
  checkmate::assertList(params)

  doc_topic <- t(lda_res$cell_topic)
  term_topic <- lda_res$topic_region

  k <- ncol(doc_topic)
  topic_names <- sprintf("topic_%02d", seq_len(k))

  rownames(doc_topic) <- doc_ids
  colnames(doc_topic) <- topic_names
  rownames(term_topic) <- term_ids
  colnames(term_topic) <- topic_names

  res <- list(
    doc_topic = doc_topic,
    term_topic = term_topic,
    bound = lda_res$bound,
    perplexity = lda_res$perplexity,
    n_iter = lda_res$n_iter,
    converged = lda_res$converged,
    doc_ids = doc_ids,
    term_ids = term_ids,
    params = params
  )

  class(res) <- "LdaResult"
  res
}

### primitives -----------------------------------------------------------------

#' @export
#'
#' @keywords internal
print.LdaResult <- function(x, ...) {
  cat("LdaResult (latent Dirichlet allocation)\n")
  cat(sprintf("  Documents:        %d\n", nrow(x$doc_topic)))
  cat(sprintf("  Terms:            %d\n", nrow(x$term_topic)))
  cat(sprintf("  Topics:           %d\n", ncol(x$doc_topic)))
  cat(sprintf("  Bound (ELBO):     %.6g\n", x$bound))
  cat(sprintf("  Perplexity:       %.6g\n", x$perplexity))
  cat(sprintf(
    "  Iterations:       %d%s\n",
    x$n_iter,
    if (x$converged) "" else " (did not converge)"
  ))
  cat("\n")
  invisible(x)
}

#' @export
#'
#' @keywords internal
dim.LdaResult <- function(x) {
  c(nrow(x$doc_topic), nrow(x$term_topic), ncol(x$doc_topic))
}

### transforms -----------------------------------------------------------------

#' @export
as.matrix.LdaResult <- function(x, which = c("doc_topic", "term_topic"), ...) {
  which <- match.arg(which)
  x[[which]]
}

### getters --------------------------------------------------------------------

#' @rdname get_params
#'
#' @export
get_params.LdaResult <- function(object, ...) {
  checkmate::assertClass(object, "LdaResult")

  object$params
}

#' @method get_params LdaResult
#'
#' @export
S7::method(get_params, S7::new_S3_class("LdaResult")) <-
  function(object, to_json = FALSE, pretty_json = FALSE) {
    get_params.LdaResult(object = object)
  }

#' Get the highest-probability terms per topic
#'
#' @description
#' The interpretability surface of a topic model: for each topic, the terms
#' carrying the most probability mass. Note that these are probabilities within
#' a topic, so they are comparable down a topic but not across topics of
#' different breadth.
#'
#' @param x `LdaResult` object.
#' @param n Integer. Number of terms to return per topic.
#'
#' @returns A data.table with `topic`, `rank`, `term` and `probability`, sorted
#' by topic and then rank.
#'
#' @export
get_top_terms <- function(x, n = 20L) {
  UseMethod("get_top_terms")
}

#' @rdname get_top_terms
#'
#' @export
get_top_terms.LdaResult <- function(x, n = 20L) {
  # checks
  checkmate::assertClass(x, "LdaResult")
  checkmate::qassert(n, "I1[1,)")

  mat <- x$term_topic
  n <- min(n, nrow(mat))

  res <- purrr::map(colnames(mat), \(topic) {
    ord <- order(mat[, topic], decreasing = TRUE)[seq_len(n)]
    data.table::data.table(
      topic = topic,
      rank = seq_len(n),
      term = rownames(mat)[ord],
      probability = mat[ord, topic]
    )
  })

  data.table::rbindlist(res)
}

## LdaKSweepResult -------------------------------------------------------------

#' Generate an LdaKSweepResult instance
#'
#' @description
#' Takes the raw Rust output of [bixverse::rs_lda_k_sweep()] and turns the
#' per-`k` metrics into a data.table, keeping the fitted models as an attribute
#' so the winning one can be pulled out without refitting.
#'
#' @param sweep_res List. The raw return of [bixverse::rs_lda_k_sweep()].
#' @param doc_ids Character vector. The document identifiers.
#' @param term_ids Character vector. The term identifiers.
#' @param params List. The parameters the sweep was run with.
#'
#' @returns An object of class `LdaKSweepResult`, which is also a data.table.
#' `combined_score` is `NA` for any `k` the coherence topic-count floor
#' excluded from selection.
#'
#' @references Arun, et al., PAKDD, 2010; Cao, et al., Neurocomputing, 2009;
#' Mimno, et al., EMNLP, 2011
#'
#' @export
#'
#' @keywords internal
new_lda_k_sweep_result <- function(sweep_res, doc_ids, term_ids, params) {
  # checks
  checkmate::assertList(sweep_res)
  checkmate::assertNames(
    names(sweep_res),
    must.include = c("k", "models", "metrics", "combined_score", "best_k")
  )
  checkmate::qassert(doc_ids, "S+")
  checkmate::qassert(term_ids, "S+")
  checkmate::assertList(params)

  models <- purrr::map(
    sweep_res$models,
    \(m) {
      new_lda_result(
        lda_res = m,
        doc_ids = doc_ids,
        term_ids = term_ids,
        params = params
      )
    }
  )
  names(models) <- as.character(sweep_res$k)

  dt <- data.table::data.table(
    k = as.integer(sweep_res$k),
    arun_2010 = purrr::map_dbl(sweep_res$metrics, "arun_2010"),
    cao_juan_2009 = purrr::map_dbl(sweep_res$metrics, "cao_juan_2009"),
    mimno_2011 = purrr::map_dbl(sweep_res$metrics, "mimno_2011"),
    bound = purrr::map_dbl(sweep_res$metrics, "bound"),
    perplexity = purrr::map_dbl(sweep_res$metrics, "perplexity"),
    combined_score = as.numeric(sweep_res$combined_score),
    converged = purrr::map_lgl(sweep_res$models, "converged")
  )

  # Rust hands back NaN for entries the coherence floor excluded. Make that an
  # honest NA rather than something that plots as a gap nobody noticed.
  excluded <- which(is.nan(dt$combined_score))
  if (length(excluded) == nrow(dt)) {
    data.table::set(dt, j = "combined_score", value = NA_real_)
    warning(sprintf(
      paste(
        "Every k in the sweep sits below the coherence topic-count floor of",
        "%d, so none of them was eligible for selection and `best_k` is not",
        "meaningful. Read the raw metrics instead, or sweep a range that",
        "reaches the floor."
      ),
      LDA_MIN_TOPICS_COH
    ))
  } else if (length(excluded) > 0L) {
    data.table::set(dt, i = excluded, j = "combined_score", value = NA_real_)
    warning(sprintf(
      paste(
        "%d of %d values of k sit below the coherence topic-count floor of %d",
        "and were excluded from the selection, so their combined_score is NA.",
        "`best_k` can therefore never be below %d, even where a smaller k has",
        "the better raw metrics. Their metrics are still reported."
      ),
      length(excluded),
      nrow(dt),
      LDA_MIN_TOPICS_COH,
      LDA_MIN_TOPICS_COH
    ))
  }

  n_unconverged <- sum(!dt$converged)
  if (n_unconverged > 0L) {
    warning(sprintf(
      paste(
        "%d of %d fits hit max_iter without the bound settling. Raise",
        "`max_iter` in params_lda() or loosen `tol` before reading too much",
        "into their metrics."
      ),
      n_unconverged,
      nrow(dt)
    ))
  }

  data.table::setattr(dt, "models", models)
  data.table::setattr(dt, "best_k", as.integer(sweep_res$best_k))
  data.table::setattr(dt, "params", params)
  data.table::setattr(
    dt,
    "class",
    c("LdaKSweepResult", "data.table", "data.frame")
  )
  dt
}

### primitives -----------------------------------------------------------------

#' @export
#'
#' @keywords internal
print.LdaKSweepResult <- function(x, ...) {
  cat("LdaKSweepResult (LDA topic count sweep)\n")
  cat(sprintf("  k range:          %d to %d\n", min(x$k), max(x$k)))
  cat(sprintf("  Best k:           %d\n", attr(x, "best_k")))
  cat(sprintf(
    "  Metrics:          arun_2010 and cao_juan_2009 lower is better,%s\n",
    " mimno_2011 higher"
  ))
  cat("\n")
  print(data.table::as.data.table(x))
  invisible(x)
}

#' Plot the LDA topic count sweep
#'
#' @description
#' The three selection metrics and the combined score against `k`. Look for the
#' `k` where the combined score peaks, then sanity check that Cao Juan has not
#' already started climbing, which is the signal that topics are duplicating.
#'
#' @param x `LdaKSweepResult` object.
#' @param ... Ignored.
#'
#' @returns A `ggplot2` object with one panel per metric.
#'
#' @export
plot.LdaKSweepResult <- function(x, ...) {
  checkmate::assertClass(x, "LdaKSweepResult")

  metrics <- c(
    arun_2010 = "Arun 2010 (lower better)",
    cao_juan_2009 = "Cao Juan 2009 (lower better)",
    mimno_2011 = "Mimno 2011 (higher better)",
    combined_score = "Combined score (higher better)"
  )

  plot_dt <- data.table::rbindlist(purrr::imap(metrics, \(label, col) {
    data.table::data.table(k = x$k, value = x[[col]], metric = label)
  }))
  plot_dt <- plot_dt[which(!is.na(plot_dt$value)), ]

  ggplot2::ggplot(
    data = plot_dt,
    mapping = ggplot2::aes(x = k, y = value)
  ) +
    ggplot2::geom_line(colour = "grey40") +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_vline(
      xintercept = attr(x, "best_k"),
      linetype = "dashed",
      colour = "firebrick"
    ) +
    ggplot2::facet_wrap(~metric, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = unique(plot_dt$k)) +
    ggplot2::labs(
      title = "LDA topic count sweep",
      subtitle = sprintf(
        "Dashed line marks the selected k = %d",
        attr(x, "best_k")
      ),
      x = "Number of topics (k)",
      y = NULL
    ) +
    ggplot2::theme_minimal()
}

### getters --------------------------------------------------------------------

#' @rdname get_params
#'
#' @export
get_params.LdaKSweepResult <- function(object, ...) {
  checkmate::assertClass(object, "LdaKSweepResult")

  attr(object, "params")
}

#' @method get_params LdaKSweepResult
#'
#' @export
S7::method(get_params, S7::new_S3_class("LdaKSweepResult")) <-
  function(object, to_json = FALSE, pretty_json = FALSE) {
    get_params.LdaKSweepResult(object = object)
  }

#' Get the selected model from an LDA topic count sweep
#'
#' @description
#' Returns the fit at `best_k`, or at a `k` you name, without refitting. The
#' sweep keeps every model it fitted.
#'
#' @details
#' `best_k` is never below five, see [bixverse::lda_k_sweep()]. Pass `k`
#' explicitly if the raw metrics point somewhere the selection could not go.
#'
#' @param x `LdaKSweepResult` object.
#' @param k Optional integer. The topic count to extract. If `NULL`, uses the
#' `best_k` the sweep selected.
#'
#' @returns An `LdaResult`.
#'
#' @export
get_best_model <- function(x, k = NULL) {
  UseMethod("get_best_model")
}

#' @rdname get_best_model
#'
#' @export
get_best_model.LdaKSweepResult <- function(x, k = NULL) {
  # checks
  checkmate::assertClass(x, "LdaKSweepResult")
  checkmate::qassert(k, c("I1", "0"))

  k <- if (is.null(k)) attr(x, "best_k") else as.integer(k)

  models <- attr(x, "models")
  if (!(as.character(k) %in% names(models))) {
    stop(sprintf(
      "No fit for k = %d in this sweep. Available: %s.",
      k,
      paste(names(models), collapse = ", ")
    ))
  }

  models[[as.character(k)]]
}
