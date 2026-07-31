use bixverse_rs::methods::nmf_bulk::{nmf_multiple_run, nmf_single_run};
use bixverse_rs::methods::nmf_hals::HalsOpts;
use bixverse_rs::prelude::*;
use extendr_api::prelude::*;

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_nmf_bulk;
    fn rs_nmf_single_bulk;
    fn rs_nmf_multi_bulk;
}

///////////////
// Functions //
///////////////

/// Run NMF (HALS) on a bulk expression matrix (single run)
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Runs a single HALS-NMF fit on a dense matrix. Expects samples x features.
/// The resulting decomposition `V ~ W H` places `W` (samples x k) as
/// sample-side factors and `H` (k x features) as feature loadings.
///
/// @param x Numerical matrix. Rows = samples, columns = features.
/// @param k Integer. Number of latent factors.
/// @param preprocessing String. One of `c("none", "sd", "sqrt_sd")`.
/// @param nmf_hals_params Named list. See [bixverse::params_nmf_hals()].
/// @param seed Integer. Random seed for the NMF initialisation.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with the following items
/// \itemize{
///   \item w - The `W` matrix of shape `n_samples x k`.
///   \item h - The `H` matrix of shape `k x n_features`.
///   \item final_loss - Final reconstruction loss.
///   \item n_iter - Number of iterations the algorithm ran for.
///   \item converged - Did the NMF algorithm converge.
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_nmf_single_bulk(
    x: RMatrix<f64>,
    k: usize,
    preprocessing: &str,
    nmf_hals_params: List,
    seed: usize,
    verbose: usize,
) -> Result<List, extendr_api::Error> {
    let x = r_matrix_to_faer(&x);
    let nmf_hals_opt: HalsOpts<f64> =
        HalsOpts::from_r_list(nmf_hals_params, seed).to_extendr()?;
    let nmf_res = nmf_single_run(x, k, preprocessing, Some(nmf_hals_opt), verbose)
        .to_extendr()?;

    Ok(list!(
        w = faer_to_r_matrix(nmf_res.w.as_ref()),
        h = faer_to_r_matrix(nmf_res.h.as_ref()),
        final_loss = nmf_res.final_loss,
        n_iter = nmf_res.n_iter,
        converged = nmf_res.converged
    ))
}

/// Run multiple NMF (HALS) restarts on a bulk expression matrix
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Runs `n_runs` HALS-NMF fits with random initialisations seeded by
/// `seed + i`. Expects samples x features.
///
/// @param x Numerical matrix. Rows = samples, columns = features.
/// @param k Integer. Number of latent factors per run.
/// @param preprocessing String. One of `c("none", "sd", "sqrt_sd")`.
/// @param nmf_hals_params Named list. See [bixverse::params_nmf_hals()].
/// @param n_runs Integer. Number of random restarts.
/// @param seed Integer. Base random seed. Run `i` uses `seed + i`.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with the following items
/// \itemize{
///   \item w_all - Column-bound `W` matrices across all runs, shape
///   `n_samples x (k * n_runs)`.
///   \item h_per_run - List of `H` matrices, each `k x n_features`.
///   \item losses - Numeric vector. Final reconstruction loss per run.
///   \item converged - Logical vector. Convergence flag per run.
///   \item best_idx - Integer. 1-indexed position of the run with the lowest
///   final loss.
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_nmf_multi_bulk(
    x: RMatrix<f64>,
    k: usize,
    preprocessing: &str,
    nmf_hals_params: List,
    n_runs: usize,
    seed: usize,
    verbose: usize,
) -> Result<List, extendr_api::Error> {
    let x = r_matrix_to_faer(&x);
    let nmf_hals_opt: HalsOpts<f64> =
        HalsOpts::from_r_list(nmf_hals_params, seed).to_extendr()?;
    let nmf_res = nmf_multiple_run(
        x,
        k,
        preprocessing,
        Some(nmf_hals_opt),
        n_runs,
        seed,
        verbose,
    )
    .to_extendr()?;

    let h_per_run: List = nmf_res
        .h_per_run
        .iter()
        .map(|h| faer_to_r_matrix(h.as_ref()))
        .collect();

    Ok(list!(
        w_all = faer_to_r_matrix(nmf_res.w_all.as_ref()),
        h_per_run = h_per_run,
        losses = nmf_res.losses,
        converged = nmf_res.converged,
        best_idx = (nmf_res.best_idx + 1) as i32
    ))
}
