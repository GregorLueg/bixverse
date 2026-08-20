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

  modules <- modules_from_loadings(gene_loadings, membership_params)

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

  modules <- modules_from_loadings(best_gene_loadings, membership_params)

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

## consensus -------------------------------------------------------------------

#' @title Run consensus NMF on a BulkCoExp
#'
#' @description
#' Runs `n_runs` HALS-NMF restarts, pools their components, drops unstable ones
#' by local density, k-means clusters the survivors into `k` groups and refits
#' the partner factor against the per-cluster median. This is cNMF: the answer
#' is the structure the restarts agree on, and the mean silhouette of those
#' clusters (`stability`) tells you how much they agreed.
#'
#' Prefer this over [bixverse::stabilised_nmf_bulk()], which just picks the
#' lowest-loss restart and tells you nothing about whether that restart is
#' reproducible.
#'
#' @details
#' Use [bixverse::nmf_k_sweep_bulk()] first if you do not already know `k`.
#'
#' The density filter is the part that bites. Components whose mean cosine
#' distance to their neighbours exceeds `density_threshold` are dropped, and if
#' that leaves fewer than `k` survivors the fit errors rather than returning
#' something degenerate. With few restarts the filter is jumpy, so either raise
#' `n_runs` or set `density_threshold = 2` to switch it off.
#'
#' @param object The class, see [bixverse::BulkCoExp()]. Ideally, you should
#' run [bixverse::preprocess_bulk_coexp()] before applying this function. The
#' NMF requires non-negative input, so the pre-processing should either skip
#' scaling or use a strictly non-negative representation.
#' @param k Integer. Number of latent factors (modules). Must be at least 2.
#' @param n_runs Integer. Number of random restarts. Must be at least 2.
#' @param preprocessing String. One of `c("none", "sd", "sqrt_sd")`. Applied
#' inside the Rust kernel to the input matrix before the HALS updates.
#' @param nmf_hals_params List. Output of [bixverse::params_nmf_hals()]. The
#' `nmf_init` field is ignored, restarts always use random initialisation.
#' @param nmf_consensus_params List. Output of
#' [bixverse::params_nmf_consensus()]. Controls the clustering of the pooled
#' components.
#' @param membership_params List. Controls how the gene loadings are turned into
#' module membership, see [bixverse::params_module_membership()].
#' @param seed Integer. Base random seed. Restart `i` uses `seed + i`, and the
#' k-means step is seeded from it too.
#' @param .verbose Boolean or integer `0L`/`1L`/`2L`. Controls verbosity.
#'
#' @return The class with `final_results` populated from the consensus fit,
#' with the clustering diagnostics under `diagnostics` and the fit parameters
#' under `params$nmf_fit`.
#'
#' @references Kotliar et al., eLife, 2019
#'
#' @export
consensus_nmf_bulk <- S7::new_generic(
  name = "consensus_nmf_bulk",
  dispatch_args = "object",
  fun = function(
    object,
    k,
    n_runs = 30L,
    preprocessing = c("none", "sd", "sqrt_sd"),
    nmf_hals_params = params_nmf_hals(),
    nmf_consensus_params = params_nmf_consensus(),
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
#' @method consensus_nmf_bulk BulkCoExp
S7::method(consensus_nmf_bulk, BulkCoExp) <- function(
  object,
  k,
  n_runs = 30L,
  preprocessing = c("none", "sd", "sqrt_sd"),
  nmf_hals_params = params_nmf_hals(),
  nmf_consensus_params = params_nmf_consensus(),
  membership_params = params_module_membership(),
  seed = 42L,
  .verbose = TRUE
) {
  preprocessing <- match.arg(preprocessing)

  # checks
  checkmate::assertClass(object, "bixverse::BulkCoExp")
  checkmate::qassert(k, "I1[2,)")
  checkmate::qassert(n_runs, "I1[2,)")
  checkmate::assertChoice(preprocessing, c("none", "sd", "sqrt_sd"))
  assertNmfHals(nmf_hals_params)
  assertNmfConsensus(nmf_consensus_params)
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
      "the data (e.g. shift, log1p) before running consensus_nmf_bulk()."
    ))
  }

  nmf_res <- .run_consensus_nmf(
    .rs_call = rs_nmf_consensus_bulk,
    nmf_consensus_params = nmf_consensus_params,
    seed = seed,
    x = target_mat,
    k = k,
    preprocessing = preprocessing,
    nmf_hals_params = nmf_hals_params,
    n_runs = n_runs,
    verbose = parse_verbosity(.verbose)
  )

  sample_activity <- nmf_res$w
  gene_loadings <- t(nmf_res$h)

  rownames(sample_activity) <- rownames(target_mat)
  rownames(gene_loadings) <- colnames(target_mat)
  colnames(sample_activity) <- colnames(gene_loadings) <-
    sprintf("comp_%02i", seq_len(k))

  modules <- modules_from_loadings(gene_loadings, membership_params)

  fit_params <- c(
    nmf_hals_params,
    list(
      k = k,
      preprocessing = preprocessing,
      seed = seed,
      membership_params = membership_params,
      nmf_consensus_params = nmf_consensus_params,
      consensus = TRUE,
      n_runs = n_runs
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
      stability = nmf_res$stability,
      rel_error = nmf_res$rel_error,
      rel_run_errors = nmf_res$rel_run_errors,
      clusters = .nmf_cluster_dt(nmf_res, k = k, n_runs = n_runs),
      cluster_sizes = data.table::data.table(
        cluster = seq_len(k),
        n = as.integer(nmf_res$cluster_sizes)
      ),
      n_dropped = nmf_res$n_dropped,
      n_empty_clusters = nmf_res$n_empty_clusters
    )
  )

  return(object)
}

## k sweep ---------------------------------------------------------------------

#' @title Sweep k for consensus NMF on a BulkCoExp
#'
#' @description
#' Runs the consensus clustering step across a range of `k` and reports
#' stability against reconstruction error, without keeping any factors. Pick
#' the `k` where stability is still high and the error curve has not yet
#' flattened out, then run [bixverse::consensus_nmf_bulk()] there.
#'
#' @details
#' This is a diagnostic, so it leaves the object alone and hands you the result
#' back directly. `plot()` on the returned object gives you the two curves.
#'
#' Cost is `length(k_range) * n_runs` full NMF fits, so keep both modest on a
#' first pass.
#'
#' @param object The class, see [bixverse::BulkCoExp()].
#' @param k_range Integer vector. The ranks to evaluate. Every entry must be at
#' least 2.
#' @param n_runs Integer. Number of random restarts per `k`. Must be at least 2.
#' @param preprocessing String. One of `c("none", "sd", "sqrt_sd")`.
#' @param nmf_hals_params List. Output of [bixverse::params_nmf_hals()].
#' @param nmf_consensus_params List. Output of
#' [bixverse::params_nmf_consensus()].
#' @param seed Integer. Base random seed.
#' @param .verbose Boolean or integer `0L`/`1L`/`2L`. Controls verbosity.
#'
#' @return An `NmfKSweepResult`, which is a data.table with one row per `k`.
#'
#' @references Kotliar et al., eLife, 2019
#'
#' @export
nmf_k_sweep_bulk <- S7::new_generic(
  name = "nmf_k_sweep_bulk",
  dispatch_args = "object",
  fun = function(
    object,
    k_range,
    n_runs = 30L,
    preprocessing = c("none", "sd", "sqrt_sd"),
    nmf_hals_params = params_nmf_hals(),
    nmf_consensus_params = params_nmf_consensus(),
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
#' @method nmf_k_sweep_bulk BulkCoExp
S7::method(nmf_k_sweep_bulk, BulkCoExp) <- function(
  object,
  k_range,
  n_runs = 30L,
  preprocessing = c("none", "sd", "sqrt_sd"),
  nmf_hals_params = params_nmf_hals(),
  nmf_consensus_params = params_nmf_consensus(),
  seed = 42L,
  .verbose = TRUE
) {
  preprocessing <- match.arg(preprocessing)

  # checks
  checkmate::assertClass(object, "bixverse::BulkCoExp")
  k_range <- .assert_nmf_k_range(k_range)
  checkmate::qassert(n_runs, "I1[2,)")
  checkmate::assertChoice(preprocessing, c("none", "sd", "sqrt_sd"))
  assertNmfHals(nmf_hals_params)
  assertNmfConsensus(nmf_consensus_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  target_mat <- .get_bulk_target_mat(object)

  if (any(target_mat < 0, na.rm = TRUE)) {
    stop(paste(
      "NMF requires non-negative input, but the target matrix contains",
      "negative values. Skip scaling in preprocess_bulk_coexp() or transform",
      "the data (e.g. shift, log1p) before running nmf_k_sweep_bulk()."
    ))
  }

  sweep_res <- rs_nmf_k_sweep_bulk(
    x = target_mat,
    k_range = k_range,
    preprocessing = preprocessing,
    nmf_hals_params = nmf_hals_params,
    nmf_consensus_params = .inject_consensus_seed(nmf_consensus_params, seed),
    n_runs = n_runs,
    seed = seed,
    verbose = parse_verbosity(.verbose)
  )

  new_nmf_k_sweep_result(
    sweep_res = sweep_res,
    source_class = "BulkCoExp",
    params = c(
      nmf_hals_params,
      list(
        k_range = k_range,
        preprocessing = preprocessing,
        seed = seed,
        nmf_consensus_params = nmf_consensus_params,
        n_runs = n_runs
      )
    )
  )
}

## helpers ---------------------------------------------------------------------

#' Validate a k range for a consensus NMF sweep
#'
#' @description
#' Checked on the R side because the disk-backed single cell sweep loads the
#' whole count matrix before Rust gets to validate anything, so a typo would
#' otherwise cost a full load.
#'
#' @param k_range Integer vector. The ranks to evaluate.
#'
#' @returns The range as an integer vector.
#'
#' @keywords internal
.assert_nmf_k_range <- function(k_range) {
  checkmate::assertIntegerish(
    k_range,
    lower = 2L,
    min.len = 1L,
    any.missing = FALSE,
    unique = TRUE
  )
  as.integer(k_range)
}

#' Seed the consensus k-means from the run seed
#'
#' @description
#' The Rust side reads a `consensus_seed` key for the k-means step.
#' [bixverse::params_nmf_consensus()] deliberately does not expose it: one call
#' having a restart seed and a separate clustering seed is a trap, where setting
#' `seed` leaves the k-means pinned at its default. So it is filled in here from
#' whatever `seed` the caller gave.
#'
#' @param nmf_consensus_params List. Output of
#' [bixverse::params_nmf_consensus()].
#' @param seed Integer. The run seed.
#'
#' @returns The parameter list with `consensus_seed` set.
#'
#' @keywords internal
.inject_consensus_seed <- function(nmf_consensus_params, seed) {
  assertNmfConsensus(nmf_consensus_params)
  checkmate::qassert(seed, "I1")
  nmf_consensus_params[["consensus_seed"]] <- as.integer(seed)
  nmf_consensus_params
}

#' Warn when clustering in cell space is about to get expensive
#'
#' @description
#' `consensus_target = "w"` clusters in sample space. On bulk that is cheap,
#' since there are tens of samples. On single cell the samples are cells, so it
#' pools a dense `(k * n_runs) x n_cells` matrix and runs an exhaustive cosine
#' search over rows of that width. Worth a heads-up before someone waits an
#' hour for it.
#'
#' @param nmf_consensus_params List. Output of
#' [bixverse::params_nmf_consensus()].
#' @param n_samples Integer. Number of cells (or meta cells) in the run.
#' @param k Integer. Components per restart.
#' @param n_runs Integer. Number of restarts.
#'
#' @returns `NULL`, invisibly. Called for the warning.
#'
#' @keywords internal
.warn_consensus_target_w <- function(
  nmf_consensus_params,
  n_samples,
  k,
  n_runs
) {
  assertNmfConsensus(nmf_consensus_params)
  checkmate::qassert(n_samples, "I1[1,)")

  cell_space_cutoff <- 20000L
  if (
    nmf_consensus_params$consensus_target == "w" &&
      n_samples > cell_space_cutoff
  ) {
    warning(sprintf(
      paste(
        "consensus_target = 'w' clusters in cell space over %d cells, which",
        "means a dense %d x %d matrix and an exhaustive cosine search across",
        "it. Use consensus_target = 'h' unless you specifically want the",
        "cell-side consensus."
      ),
      n_samples,
      k * n_runs,
      n_samples
    ))
  }

  invisible(NULL)
}

#' Run a consensus NMF binding with an actionable error
#'
#' @description
#' The consensus step can only fail once the restarts are done: the density
#' filter may leave fewer than `k` components, or a cluster may end up empty.
#' Both are real signals that the factorisation is not reproducible at this `k`,
#' so neither is retried or degraded into a partial answer. The Rust message
#' says what happened; this appends what to do about it.
#'
#' @param .rs_call Function. The `rs_nmf_consensus_*` binding to call.
#' @param nmf_consensus_params List. Output of
#' [bixverse::params_nmf_consensus()]. The consensus seed is injected here.
#' @param seed Integer. The run seed.
#' @param ... Remaining arguments forwarded to `.rs_call`.
#'
#' @returns The raw list the binding returned.
#'
#' @keywords internal
.run_consensus_nmf <- function(.rs_call, nmf_consensus_params, seed, ...) {
  checkmate::assertFunction(.rs_call)
  assertNmfConsensus(nmf_consensus_params)
  checkmate::qassert(seed, "I1")

  tryCatch(
    .rs_call(
      nmf_consensus_params = .inject_consensus_seed(nmf_consensus_params, seed),
      seed = seed,
      ...
    ),
    error = function(e) {
      stop(paste(
        conditionMessage(e),
        "\nOptions: raise `density_threshold` in params_nmf_consensus()",
        "(2 switches the filter off entirely), increase `n_runs`, or run the",
        "k sweep first to find a k the restarts actually agree on."
      ))
    }
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

#' @title Get the multi-run NMF diagnostics
#'
#' @description
#' Getter function for the diagnostics of a multi-run NMF fit. For a
#' [bixverse::stabilised_nmf_bulk()] fit that is the per-run losses,
#' convergence flags and the best-run index. For a
#' [bixverse::consensus_nmf_bulk()] fit it is the cluster stability, the
#' relative errors and the per-component clustering table. Returns `NULL` for a
#' single-run fit.
#'
#' @param object The class, see [bixverse::BulkCoExp()].
#'
#' @return For a stabilised fit, a list with `losses`, `converged` and
#' `best_idx`. For a consensus fit, a list with `stability`, `rel_error`,
#' `rel_run_errors`, `clusters`, `cluster_sizes`, `n_dropped` and
#' `n_empty_clusters`. `NULL` if neither was run.
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

  no_diagnostics <- function() {
    warning(paste(
      "No multi-run NMF diagnostics found. Did you run",
      "stabilised_nmf_bulk() or consensus_nmf_bulk()? Returning NULL."
    ))
    NULL
  }

  final_results <- S7::prop(object, "final_results")
  if (!inherits(final_results, "BulkModuleResult")) {
    return(no_diagnostics())
  }
  diagnostics <- get_diagnostics(final_results)

  if (!is.null(diagnostics[["stability"]])) {
    return(list(
      stability = diagnostics[["stability"]],
      rel_error = diagnostics[["rel_error"]],
      rel_run_errors = diagnostics[["rel_run_errors"]],
      clusters = diagnostics[["clusters"]],
      cluster_sizes = diagnostics[["cluster_sizes"]],
      n_dropped = diagnostics[["n_dropped"]],
      n_empty_clusters = diagnostics[["n_empty_clusters"]]
    ))
  }

  if (is.null(diagnostics[["losses"]])) {
    return(no_diagnostics())
  }
  list(
    losses = diagnostics[["losses"]],
    converged = diagnostics[["converged"]],
    best_idx = diagnostics[["best_idx"]]
  )
}
