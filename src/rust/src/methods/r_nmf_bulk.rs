use bixverse_rs::methods::nmf_bulk::{
    nmf_consensus_run, nmf_k_sweep_run, nmf_multiple_run, nmf_single_run,
};
use bixverse_rs::methods::nmf_hals::consensus::ConsensusParams;
use bixverse_rs::methods::nmf_hals::HalsOpts;
use bixverse_rs::prelude::*;
use extendr_api::prelude::*;

use crate::methods::nmf_utils::{consensus_res_to_r_list, k_sweep_to_r_list};

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_nmf_bulk;
    fn rs_nmf_single_bulk;
    fn rs_nmf_multi_bulk;
    fn rs_nmf_consensus_bulk;
    fn rs_nmf_k_sweep_bulk;
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

/// Run consensus NMF on a bulk expression matrix
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Runs `n_runs` HALS-NMF restarts, pools their components, drops unstable
/// ones by local density, k-means clusters the survivors and refits the
/// partner factor against the per-cluster median. Expects samples x features.
///
/// @param x Numerical matrix. Rows = samples, columns = features.
/// @param k Integer. Number of latent factors. Must be at least 2.
/// @param preprocessing String. One of `c("none", "sd", "sqrt_sd")`.
/// @param nmf_hals_params Named list. See [bixverse::params_nmf_hals()]. The
/// `nmf_init` field is ignored, restarts always use random initialisation.
/// @param nmf_consensus_params Named list. See
/// [bixverse::params_nmf_consensus()].
/// @param n_runs Integer. Number of restarts. Must be at least 2.
/// @param seed Integer. Base random seed. Restart `i` uses `seed + i`.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with the following items
/// \itemize{
///   \item w - The `W` matrix of shape `n_samples x k`.
///   \item h - The `H` matrix of shape `k x n_features`.
///   \item rel_error - Reconstruction error relative to the squared Frobenius
///   norm of the input. Not comparable with the absolute `final_loss` the
///   single-run version returns.
///   \item rel_run_errors - The same, per restart.
///   \item labels - Integer vector of length `k * n_runs`. Cluster each pooled
///   component landed in, `NA` if it was dropped.
///   \item local_density - Mean cosine distance to the nearest neighbours per
///   pooled component. Plot this to pick a `density_threshold`.
///   \item kept - 1-indexed positions of the surviving pooled components.
///   \item silhouette - Silhouette per survivor, aligned with `kept`.
///   \item stability - Mean silhouette over the survivors.
///   \item cluster_sizes - Number of survivors per cluster.
///   \item n_dropped - Number of pooled components removed.
///   \item n_empty_clusters - Number of clusters left with no members.
/// }
///
/// @references Kotliar et al., eLife, 2019
///
/// @export
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_nmf_consensus_bulk(
    x: RMatrix<f64>,
    k: usize,
    preprocessing: &str,
    nmf_hals_params: List,
    nmf_consensus_params: List,
    n_runs: usize,
    seed: usize,
    verbose: usize,
) -> Result<List, extendr_api::Error> {
    let x = r_matrix_to_faer(&x);
    let nmf_hals_opt: HalsOpts<f64> =
        HalsOpts::from_r_list(nmf_hals_params, seed).to_extendr()?;
    let consensus_opt: ConsensusParams<f64> =
        ConsensusParams::from_r_list(nmf_consensus_params)?;
    let nmf_res = nmf_consensus_run(
        x,
        k,
        preprocessing,
        Some(nmf_hals_opt),
        Some(consensus_opt),
        n_runs,
        seed,
        verbose,
    )
    .to_extendr()?;

    Ok(consensus_res_to_r_list(&nmf_res))
}

/// Sweep k and report consensus stability against reconstruction error (bulk)
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Returns diagnostics only, no factors, so a wide `k_range` stays cheap.
/// Pick the k where stability is high and the error curve has not yet
/// flattened, then call [rs_nmf_consensus_bulk()] there.
///
/// @param x Numerical matrix. Rows = samples, columns = features.
/// @param k_range Integer vector. Ranks to evaluate, every entry at least 2.
/// @param preprocessing String. One of `c("none", "sd", "sqrt_sd")`.
/// @param nmf_hals_params Named list. See [bixverse::params_nmf_hals()].
/// @param nmf_consensus_params Named list. See
/// [bixverse::params_nmf_consensus()].
/// @param n_runs Integer. Number of restarts per k. Must be at least 2.
/// @param seed Integer. Base random seed.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list of equal-length vectors, one element per swept k
/// \itemize{
///   \item k - The rank.
///   \item stability - Mean silhouette of the consensus clusters. `NaN` where
///   the consensus step failed.
///   \item best_error - Lowest restart error, relative to the squared
///   Frobenius norm of the input.
///   \item median_error - Median restart error, same scale.
///   \item consensus_failed - Did the density filter leave fewer than `k`
///   components.
///   \item n_dropped - Number of pooled components removed.
///   \item n_empty_clusters - Number of clusters left with no members.
///   \item n_converged - Restarts that met the HALS tolerance.
/// }
///
/// @references Kotliar et al., eLife, 2019
///
/// @export
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_nmf_k_sweep_bulk(
    x: RMatrix<f64>,
    k_range: &[i32],
    preprocessing: &str,
    nmf_hals_params: List,
    nmf_consensus_params: List,
    n_runs: usize,
    seed: usize,
    verbose: usize,
) -> Result<List, extendr_api::Error> {
    let x = r_matrix_to_faer(&x);
    let k_range = k_range.r_int_convert();
    let nmf_hals_opt: HalsOpts<f64> =
        HalsOpts::from_r_list(nmf_hals_params, seed).to_extendr()?;
    let consensus_opt: ConsensusParams<f64> =
        ConsensusParams::from_r_list(nmf_consensus_params)?;
    let sweep_res = nmf_k_sweep_run(
        x,
        &k_range,
        preprocessing,
        Some(nmf_hals_opt),
        Some(consensus_opt),
        n_runs,
        seed,
        verbose,
    )
    .to_extendr()?;

    Ok(k_sweep_to_r_list(&sweep_res))
}
