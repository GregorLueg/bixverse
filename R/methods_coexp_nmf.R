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
#' Module membership is derived by keeping the upper tail of each component's
#' gene loadings, see [bixverse::params_module_membership()]. A gene can belong
#' to several modules, and a gene in no tail belongs to none.
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
#' @param membership_params List. Controls how the gene loadings are turned into
#' module membership, see [bixverse::params_module_membership()]. Membership is
#' not exclusive: a gene loading strongly on several components appears in
#' several modules, and a gene in no tail appears in none.
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
    membership_params = params_module_membership(),
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
  membership_params = params_module_membership(),
  seed = 42L,
  .verbose = TRUE
) {
  preprocessing <- match.arg(preprocessing)

  # checks
  checkmate::assertClass(object, "bixverse::BulkCoExp")
  checkmate::qassert(k, "I1[1,)")
  checkmate::assertChoice(preprocessing, c("none", "sd", "sqrt_sd"))
  assertNmfHals(nmf_hals_params)
  assertModuleMembershipParams(membership_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  # early return
  detection_method <- .assert_bulk_detection_method(
    object,
    "nmf-based",
    "nmf-based",
    allow_unset = TRUE
  )
  if (is.null(detection_method)) {
    return(object)
  }

  # target matrix (samples x features)
  target_mat <- .get_bulk_target_mat(object)

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

  modules <- .nmf_modules_from_w(gene_loadings, membership_params)

  fit_params <- c(
    nmf_hals_params,
    list(
      k = k,
      preprocessing = preprocessing,
      seed = seed,
      membership_params = membership_params,
      stabilised = FALSE,
      final_loss = nmf_res$final_loss,
      n_iter = nmf_res$n_iter,
      converged = nmf_res$converged
    )
  )

  S7::prop(object, "params")[["nmf_fit"]] <- fit_params
  S7::prop(object, "params")[["detection_method"]] <- "nmf-based"

  S7::prop(object, "final_results") <- new_bulk_module_result(
    modules = modules,
    factors = list(
      gene_loadings = gene_loadings,
      sample_activity = sample_activity
    ),
    method = "nmf-based",
    params = S7::prop(object, "params"),
    diagnostics = list(
      final_loss = nmf_res$final_loss,
      n_iter = nmf_res$n_iter,
      converged = nmf_res$converged
    )
  )

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
#' @param membership_params List. Controls how the gene loadings are turned into
#' module membership, see [bixverse::params_module_membership()].
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
    membership_params = params_module_membership(),
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
  membership_params = params_module_membership(),
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
  assertModuleMembershipParams(membership_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  # early return
  detection_method <- .assert_bulk_detection_method(
    object,
    "nmf-based",
    "nmf-based",
    allow_unset = TRUE
  )
  if (is.null(detection_method)) {
    return(object)
  }

  # target matrix (samples x features)
  target_mat <- .get_bulk_target_mat(object)

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

  modules <- .nmf_modules_from_w(best_gene_loadings, membership_params)

  fit_params <- c(
    nmf_hals_params,
    list(
      k = k,
      preprocessing = preprocessing,
      seed = seed,
      membership_params = membership_params,
      stabilised = TRUE,
      n_runs = n_runs,
      best_idx = best_idx,
      final_loss = nmf_res$losses[best_idx]
    )
  )

  S7::prop(object, "params")[["nmf_fit"]] <- fit_params
  S7::prop(object, "params")[["detection_method"]] <- "nmf-based"

  S7::prop(object, "final_results") <- new_bulk_module_result(
    modules = modules,
    factors = list(
      gene_loadings = best_gene_loadings,
      sample_activity = best_sample_activity
    ),
    method = "nmf-based",
    params = S7::prop(object, "params"),
    diagnostics = list(
      losses = nmf_res$losses,
      converged = nmf_res$converged,
      best_idx = best_idx,
      w_all_runs = nmf_res$w_all,
      h_per_run = nmf_res$h_per_run
    )
  )

  return(object)
}

## helpers ---------------------------------------------------------------------

#' Derive gene-to-module assignments from an NMF loadings matrix
#'
#' @description
#' Each gene is assigned to the component whose absolute loading is highest.
#' Thin wrapper over [bixverse::.modules_from_loadings()]. NMF loadings are
#' non-negative, so the `"auto"` tail setting resolves to an upper-tail test.
#'
#' @param w Numeric matrix. Gene loadings (features x k).
#' @param membership_params List. See
#' [bixverse::params_module_membership()].
#'
#' @returns A data.table with one row per surviving (gene, component) pair.
#'
#' @keywords internal
.nmf_modules_from_w <- function(
  w,
  membership_params = params_module_membership()
) {
  .modules_from_loadings(w, membership_params)
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
  if (!inherits(final_results, "BulkModuleResult")) {
    warning(paste(
      "No NMF gene loadings found. Did you run nmf_bulk() or",
      "stabilised_nmf_bulk()? Returning NULL."
    ))
    return(NULL)
  }
  get_factors(final_results, which = "gene_loadings")
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
  if (!inherits(final_results, "BulkModuleResult")) {
    warning(paste(
      "No NMF sample activity found. Did you run nmf_bulk() or",
      "stabilised_nmf_bulk()? Returning NULL."
    ))
    return(NULL)
  }
  get_factors(final_results, which = "sample_activity")
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
  if (!inherits(final_results, "BulkModuleResult")) {
    warning(paste(
      "No NMF modules found. Did you run nmf_bulk() or",
      "stabilised_nmf_bulk()? Returning NULL."
    ))
    return(NULL)
  }
  get_modules(final_results)
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
  if (!inherits(final_results, "BulkModuleResult")) {
    warning(paste(
      "No stabilised NMF diagnostics found. Did you run",
      "stabilised_nmf_bulk()? Returning NULL."
    ))
    return(NULL)
  }
  diagnostics <- get_diagnostics(final_results)
  if (is.null(diagnostics[["losses"]])) {
    warning(paste(
      "No stabilised NMF diagnostics found. Did you run",
      "stabilised_nmf_bulk()? Returning NULL."
    ))
    return(NULL)
  }
  list(
    losses = diagnostics[["losses"]],
    converged = diagnostics[["converged"]],
    best_idx = diagnostics[["best_idx"]]
  )
}
