# methods - NMF ----------------------------------------------------------------

## single run ------------------------------------------------------------------

#' @title Run non-negative matrix factorisation on a BulkCoExp
#'
#' @description
#' Fits a single HALS-NMF `V ~ W H` to the bulk expression matrix with a fixed
#' number of components `k`. `V` is samples x features (matches the layout of
#' `raw_data` / `processed_data` on [bixverse::BulkCoExp()]). The result is
#' rearranged so that:
#' \itemize{
#'  \item `gene_loadings` (features x k) captures per-module gene contributions.
#'  \item `sample_activity` (samples x k) captures per-module sample activity.
#' }
#' Module membership is derived per gene by `which.max` over `abs(gene_loadings)`
#' with a minimum-loading threshold.
#'
#' @param object The class, see [bixverse::BulkCoExp()]. Ideally, you should
#' run [bixverse::preprocess_bulk_coexp()] before applying this function. The
#' NMF requires non-negative input, so the pre-processing should either skip
#' scaling or use a strictly non-negative representation.
#' @param k Integer. Number of latent factors (modules).
#' @param preprocessing String. One of `c("none", "sd", "sqrt_sd")`. Applied
#' inside the Rust kernel to the input matrix before the HALS updates.
#' @param nmf_hals_params List. Output of [bixverse::params_nmf_hals()]:
#' \itemize{
#'  \item max_iter - Integer. Maximum number of HALS iterations.
#'  \item tol - Float. Convergence tolerance.
#'  \item eps - Float. Numerical floor.
#'  \item check_every - Integer. Convergence check interval.
#'  \item nmf_init - String. One of `c("nndsvd", "svd", "random")`.
#' }
#' @param min_loading Float in `[0, 1]`. Minimum fraction of the top loading a
#' gene must reach to be assigned to any module. Genes falling below the
#' threshold are labelled `"unassigned"`. Defaults to `0`.
#' @param seed Integer. Random seed for the NMF initialisation. Defaults to
#' `42L`.
#' @param .verbose Boolean or integer `0L`/`1L`/`2L`. Controls verbosity.
#'
#' @return The class with `final_results` populated (see description) and the
#' fit parameters stored under `params$nmf_fit`, plus
#' `params$detection_method = "nmf-based"`.
#'
#' @references Cichocki & Phan, IEICE Trans., 2009.
#'
#' @export
nmf_bulk <- S7::new_generic(
  name = "nmf_bulk",
  dispatch_args = "object",
  fun = function(
    object,
    k,
    preprocessing = c("none", "sd", "sqrt_sd"),
    nmf_hals_params = params_nmf_hals(),
    min_loading = 0,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @import data.table
#'
#' @method nmf_bulk BulkCoExp
S7::method(nmf_bulk, BulkCoExp) <- function(
  object,
  k,
  preprocessing = c("none", "sd", "sqrt_sd"),
  nmf_hals_params = params_nmf_hals(),
  min_loading = 0,
  seed = 42L,
  .verbose = TRUE
) {
  preprocessing <- match.arg(preprocessing)

  # checks
  checkmate::assertClass(object, "bixverse::BulkCoExp")
  checkmate::qassert(k, "I1[1,)")
  checkmate::assertChoice(preprocessing, c("none", "sd", "sqrt_sd"))
  assertNmfHals(nmf_hals_params)
  checkmate::qassert(min_loading, "N1[0, 1]")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  # target matrix (samples x features)
  if (purrr::is_empty(S7::prop(object, "processed_data")[["processed_data"]])) {
    warning("No pre-processed data found. Defaulting to the raw data.")
    target_mat <- S7::prop(object, "raw_data")
  } else {
    target_mat <- S7::prop(object, "processed_data")[["processed_data"]]
  }

  if (any(target_mat < 0, na.rm = TRUE)) {
    stop(paste(
      "NMF requires non-negative input, but the target matrix contains",
      "negative values. Skip scaling in preprocess_bulk_coexp() or transform",
      "the data (e.g. shift, log1p) before running nmf_bulk()."
    ))
  }

  verbose_int <- as.integer(as.logical(.verbose))
  if (is.numeric(.verbose)) {
    verbose_int <- as.integer(.verbose)
  }

  nmf_res <- rs_nmf_single_bulk(
    x = target_mat,
    k = k,
    preprocessing = preprocessing,
    nmf_hals_params = nmf_hals_params,
    seed = seed,
    verbose = verbose_int
  )

  # w is samples x k, h is k x features (matches Rust output)
  sample_activity <- nmf_res$w
  gene_loadings <- t(nmf_res$h)

  rownames(sample_activity) <- rownames(target_mat)
  rownames(gene_loadings) <- colnames(target_mat)
  colnames(sample_activity) <- colnames(gene_loadings) <-
    sprintf("comp_%02i", seq_len(k))

  modules <- .nmf_modules_from_w(gene_loadings, min_loading)

  fit_params <- c(
    nmf_hals_params,
    list(
      k = k,
      preprocessing = preprocessing,
      seed = seed,
      min_loading = min_loading,
      stabilised = FALSE,
      final_loss = nmf_res$final_loss,
      n_iter = nmf_res$n_iter,
      converged = nmf_res$converged
    )
  )

  results <- list(
    modules = modules,
    gene_loadings = gene_loadings,
    sample_activity = sample_activity,
    final_loss = nmf_res$final_loss,
    n_iter = nmf_res$n_iter,
    converged = nmf_res$converged
  )

  S7::prop(object, "params")[["nmf_fit"]] <- fit_params
  S7::prop(object, "params")[["detection_method"]] <- "nmf-based"
  S7::prop(object, "final_results") <- results

  return(object)
}

## stabilised ------------------------------------------------------------------

#' @title Run stabilised (multi-restart) NMF on a BulkCoExp
#'
#' @description
#' Runs `n_runs` HALS-NMF fits with random initialisations and returns the run
#' with the lowest final reconstruction loss as the primary result, along with
#' the per-run losses and convergence flags. Useful when the objective is
#' known to have multiple local minima. Uses `f64` precision.
#'
#' @param object The class, see [bixverse::BulkCoExp()].
#' @param k Integer. Number of latent factors (modules) per run.
#' @param n_runs Integer. Number of random restarts.
#' @param preprocessing String. One of `c("none", "sd", "sqrt_sd")`.
#' @param nmf_hals_params List. Output of [bixverse::params_nmf_hals()]. Note
#' that `nmf_init` is ignored for the stabilised variant, which always uses
#' random initialisation seeded by `seed + i`.
#' @param min_loading Float in `[0, 1]`. See [bixverse::nmf_bulk()].
#' @param seed Integer. Base random seed.
#' @param .verbose Boolean or integer `0L`/`1L`/`2L`. Controls verbosity.
#'
#' @return The class with `final_results` populated from the best run (lowest
#' final loss) plus stability diagnostics (`losses`, `converged`, `best_idx`,
#' `w_all_runs`, `h_per_run`).
#'
#' @export
stabilised_nmf_bulk <- S7::new_generic(
  name = "stabilised_nmf_bulk",
  dispatch_args = "object",
  fun = function(
    object,
    k,
    n_runs = 30L,
    preprocessing = c("none", "sd", "sqrt_sd"),
    nmf_hals_params = params_nmf_hals(),
    min_loading = 0,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @import data.table
#'
#' @method stabilised_nmf_bulk BulkCoExp
S7::method(stabilised_nmf_bulk, BulkCoExp) <- function(
  object,
  k,
  n_runs = 30L,
  preprocessing = c("none", "sd", "sqrt_sd"),
  nmf_hals_params = params_nmf_hals(),
  min_loading = 0,
  seed = 42L,
  .verbose = TRUE
) {
  preprocessing <- match.arg(preprocessing)

  # checks
  checkmate::assertClass(object, "bixverse::BulkCoExp")
  checkmate::qassert(k, "I1[1,)")
  checkmate::qassert(n_runs, "I1[1,)")
  checkmate::assertChoice(preprocessing, c("none", "sd", "sqrt_sd"))
  assertNmfHals(nmf_hals_params)
  checkmate::qassert(min_loading, "N1[0, 1]")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  # target matrix (samples x features)
  if (purrr::is_empty(S7::prop(object, "processed_data")[["processed_data"]])) {
    warning("No pre-processed data found. Defaulting to the raw data.")
    target_mat <- S7::prop(object, "raw_data")
  } else {
    target_mat <- S7::prop(object, "processed_data")[["processed_data"]]
  }

  if (any(target_mat < 0, na.rm = TRUE)) {
    stop(paste(
      "NMF requires non-negative input, but the target matrix contains",
      "negative values. Skip scaling in preprocess_bulk_coexp() or transform",
      "the data (e.g. shift, log1p) before running stabilised_nmf_bulk()."
    ))
  }

  verbose_int <- as.integer(as.logical(.verbose))
  if (is.numeric(.verbose)) {
    verbose_int <- as.integer(.verbose)
  }

  nmf_res <- rs_nmf_multi_bulk(
    x = target_mat,
    k = k,
    preprocessing = preprocessing,
    nmf_hals_params = nmf_hals_params,
    n_runs = n_runs,
    seed = seed,
    verbose = verbose_int
  )

  best_idx <- nmf_res$best_idx
  # w_all is (samples x (k * n_runs)); pull the columns for best_idx
  block_start <- (best_idx - 1L) * k + 1L
  block_end <- best_idx * k
  best_sample_activity <- nmf_res$w_all[, block_start:block_end, drop = FALSE]
  best_h <- nmf_res$h_per_run[[best_idx]]
  best_gene_loadings <- t(best_h)

  rownames(best_sample_activity) <- rownames(target_mat)
  rownames(best_gene_loadings) <- colnames(target_mat)
  colnames(best_sample_activity) <- colnames(best_gene_loadings) <-
    sprintf("comp_%02i", seq_len(k))

  modules <- .nmf_modules_from_w(best_gene_loadings, min_loading)

  fit_params <- c(
    nmf_hals_params,
    list(
      k = k,
      preprocessing = preprocessing,
      seed = seed,
      min_loading = min_loading,
      stabilised = TRUE,
      n_runs = n_runs,
      best_idx = best_idx,
      final_loss = nmf_res$losses[best_idx]
    )
  )

  results <- list(
    modules = modules,
    gene_loadings = best_gene_loadings,
    sample_activity = best_sample_activity,
    losses = nmf_res$losses,
    converged = nmf_res$converged,
    best_idx = best_idx,
    w_all_runs = nmf_res$w_all,
    h_per_run = nmf_res$h_per_run
  )

  S7::prop(object, "params")[["nmf_fit"]] <- fit_params
  S7::prop(object, "params")[["detection_method"]] <- "nmf-based"
  S7::prop(object, "final_results") <- results

  return(object)
}

## helpers ---------------------------------------------------------------------

#' Derive gene-to-module assignments from an NMF loadings matrix
#'
#' @description
#' Each gene is assigned to the component whose absolute loading is highest.
#' Genes whose top loading is below `min_loading * max(|W|)` are labelled as
#' `"unassigned"` and excluded from the module table.
#'
#' @param w Numeric matrix. Gene loadings (features x k).
#' @param min_loading Numeric. Minimum fraction of the global max loading a
#' gene must reach to be assigned to a module.
#'
#' @returns A data.table with columns `gene`, `module_id`, `loading`, `sign`.
#'
#' @keywords internal
.nmf_modules_from_w <- function(w, min_loading = 0) {
  checkmate::assertMatrix(w, mode = "numeric", row.names = "named")
  checkmate::qassert(min_loading, "N1[0, 1]")

  abs_w <- abs(w)
  top_idx <- apply(abs_w, 1L, which.max)
  top_val <- abs_w[cbind(seq_len(nrow(w)), top_idx)]
  sign_val <- sign(w[cbind(seq_len(nrow(w)), top_idx)])

  cutoff <- min_loading * max(abs_w)
  keep <- top_val >= cutoff

  data.table::data.table(
    gene = rownames(w)[keep],
    module_id = colnames(w)[top_idx[keep]],
    loading = top_val[keep],
    sign = sign_val[keep]
  )
}

## getters ---------------------------------------------------------------------

#' @title Get the NMF gene loadings
#'
#' @description
#' Getter function to extract the gene loadings matrix (features x k) from a
#' bulk NMF fit stored in [bixverse::BulkCoExp()]. If no NMF fit is present,
#' returns `NULL` with a warning.
#'
#' @param object The class, see [bixverse::BulkCoExp()].
#'
#' @return A features x k numeric matrix (if found) or `NULL`.
#'
#' @export
get_nmf_gene_loadings <- S7::new_generic(
  name = "get_nmf_gene_loadings",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' @export
#' @method get_nmf_gene_loadings BulkCoExp
S7::method(get_nmf_gene_loadings, BulkCoExp) <- function(object) {
  checkmate::assertClass(object, "bixverse::BulkCoExp")
  final_results <- S7::prop(object, "final_results")
  loadings <- if (is.list(final_results)) final_results[["gene_loadings"]]
  if (is.null(loadings)) {
    warning(paste(
      "No NMF gene loadings found. Did you run nmf_bulk() or",
      "stabilised_nmf_bulk()? Returning NULL."
    ))
  }
  return(loadings)
}

#' @title Get the NMF sample activity
#'
#' @description
#' Getter function to extract the sample activity matrix (samples x k) from a
#' bulk NMF fit. If no NMF fit is present, returns `NULL` with a warning.
#'
#' @param object The class, see [bixverse::BulkCoExp()].
#'
#' @return A samples x k numeric matrix (if found) or `NULL`.
#'
#' @export
get_nmf_sample_activity <- S7::new_generic(
  name = "get_nmf_sample_activity",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' @export
#' @method get_nmf_sample_activity BulkCoExp
S7::method(get_nmf_sample_activity, BulkCoExp) <- function(object) {
  checkmate::assertClass(object, "bixverse::BulkCoExp")
  final_results <- S7::prop(object, "final_results")
  activity <- if (is.list(final_results)) final_results[["sample_activity"]]
  if (is.null(activity)) {
    warning(paste(
      "No NMF sample activity found. Did you run nmf_bulk() or",
      "stabilised_nmf_bulk()? Returning NULL."
    ))
  }
  return(activity)
}

#' @title Get the NMF module membership data.table
#'
#' @description
#' Getter function to extract the gene-to-module data.table from a bulk NMF
#' fit. Each row is one gene assigned to its top-loading module.
#'
#' @param object The class, see [bixverse::BulkCoExp()].
#'
#' @return A data.table with `gene`, `module_id`, `loading`, `sign` columns
#' (if found) or `NULL`.
#'
#' @export
get_nmf_modules <- S7::new_generic(
  name = "get_nmf_modules",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' @export
#' @method get_nmf_modules BulkCoExp
S7::method(get_nmf_modules, BulkCoExp) <- function(object) {
  checkmate::assertClass(object, "bixverse::BulkCoExp")
  final_results <- S7::prop(object, "final_results")
  modules <- if (is.list(final_results)) final_results[["modules"]]
  if (is.null(modules)) {
    warning(paste(
      "No NMF modules found. Did you run nmf_bulk() or",
      "stabilised_nmf_bulk()? Returning NULL."
    ))
  }
  return(modules)
}

#' @title Get the stabilised NMF diagnostics
#'
#' @description
#' Getter function to extract per-run losses, convergence flags, and the
#' best-run index from a `stabilised_nmf_bulk()` fit. Returns `NULL` for a
#' single-run fit.
#'
#' @param object The class, see [bixverse::BulkCoExp()].
#'
#' @return A list with `losses`, `converged`, `best_idx` (if found) or `NULL`.
#'
#' @export
get_nmf_stability <- S7::new_generic(
  name = "get_nmf_stability",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' @export
#' @method get_nmf_stability BulkCoExp
S7::method(get_nmf_stability, BulkCoExp) <- function(object) {
  checkmate::assertClass(object, "bixverse::BulkCoExp")
  final_results <- S7::prop(object, "final_results")
  if (!is.list(final_results) || is.null(final_results[["losses"]])) {
    warning(paste(
      "No stabilised NMF diagnostics found. Did you run",
      "stabilised_nmf_bulk()? Returning NULL."
    ))
    return(NULL)
  }
  return(list(
    losses = final_results[["losses"]],
    converged = final_results[["converged"]],
    best_idx = final_results[["best_idx"]]
  ))
}
