//! Shared conversions for the consensus NMF bindings.
//!
//! The bulk, single cell and meta cell entry points differ only in how they get
//! at the data and in the float they run at. What comes back is the same shape,
//! so the R list is built once here.
//!
//! Everything is flattened deliberately: the R constructors reshape, name and
//! bind the pieces together, the same way `new_nmf_result()` does for the
//! single-run fits.

use bixverse_rs::methods::nmf_hals::consensus::{ConsensusNmfResult, KSweepEntry};
use bixverse_rs::prelude::*;
use extendr_api::prelude::*;

/////////////////
// Conversions //
/////////////////

/// Convert a consensus NMF fit into a flat R list
///
/// Cluster labels and the pooled indices of the survivors are shifted to
/// 1-based on the way out. A dropped component has no label, which becomes
/// `NA_integer_` rather than a sentinel the R side would have to know about.
///
/// ### Params
///
/// * `res` - The fit to convert.
///
/// ### Returns
///
/// A list with `w`, `h`, `rel_error`, `rel_run_errors`, `labels`,
/// `local_density`, `kept`, `silhouette`, `stability`, `cluster_sizes`,
/// `n_dropped` and `n_empty_clusters`. `rel_error` and `rel_run_errors` are
/// relative to the squared Frobenius norm of the input, unlike the absolute
/// losses the single-run and stabilised bindings return.
pub fn consensus_res_to_r_list<F>(res: &ConsensusNmfResult<F>) -> List
where
    F: BixverseFloat + FaerRType,
{
    let labels: Vec<Rint> = res
        .clusters
        .labels
        .iter()
        .map(|label| match label {
            Some(cluster) => Rint::from((*cluster + 1) as i32),
            None => Rint::na(),
        })
        .collect();

    let kept: Vec<i32> = res
        .clusters
        .kept
        .iter()
        .map(|&idx| (idx + 1) as i32)
        .collect();

    let cluster_sizes: Vec<i32> = res.clusters.sizes.iter().map(|&s| s as i32).collect();

    list!(
        w = faer_to_r_matrix(res.w.as_ref()),
        h = faer_to_r_matrix(res.h.as_ref()),
        rel_error = to_f64(res.error),
        rel_run_errors = floats_to_f64(&res.run_errors),
        labels = labels,
        local_density = floats_to_f64(&res.clusters.local_density),
        kept = kept,
        silhouette = floats_to_f64(&res.clusters.silhouette),
        stability = to_f64(res.clusters.stability),
        cluster_sizes = cluster_sizes,
        n_dropped = res.clusters.n_dropped as i32,
        n_empty_clusters = res.clusters.n_empty_clusters as i32
    )
}

/// Convert a k sweep into a columnar R list
///
/// One element per field rather than one per k, so the R side can hand the
/// whole thing to `data.table::as.data.table()` without a transpose.
///
/// ### Params
///
/// * `entries` - One entry per swept k, in the order they were requested.
///
/// ### Returns
///
/// A list with `k`, `stability`, `best_error`, `median_error`,
/// `consensus_failed`, `n_dropped`, `n_empty_clusters` and `n_converged`, each
/// a vector as long as `entries`. `stability` is `NaN` where the consensus step
/// failed; the R side turns that into `NA`.
pub fn k_sweep_to_r_list<F>(entries: &[KSweepEntry<F>]) -> List
where
    F: BixverseFloat,
{
    let k: Vec<i32> = entries.iter().map(|e| e.k as i32).collect();
    let stability: Vec<f64> = entries.iter().map(|e| to_f64(e.stability)).collect();
    let best_error: Vec<f64> = entries.iter().map(|e| to_f64(e.best_error)).collect();
    let median_error: Vec<f64> = entries.iter().map(|e| to_f64(e.median_error)).collect();
    let consensus_failed: Vec<bool> = entries.iter().map(|e| e.consensus_failed).collect();
    let n_dropped: Vec<i32> = entries.iter().map(|e| e.n_dropped as i32).collect();
    let n_empty_clusters: Vec<i32> = entries.iter().map(|e| e.n_empty_clusters as i32).collect();
    let n_converged: Vec<i32> = entries.iter().map(|e| e.n_converged as i32).collect();

    list!(
        k = k,
        stability = stability,
        best_error = best_error,
        median_error = median_error,
        consensus_failed = consensus_failed,
        n_dropped = n_dropped,
        n_empty_clusters = n_empty_clusters,
        n_converged = n_converged
    )
}

/////////////
// Helpers //
/////////////

/// Widen a single float to `f64` for the trip across the boundary.
///
/// `NaN` survives the conversion, which the k sweep relies on to signal a
/// failed consensus step.
///
/// ### Params
///
/// * `x` - The value to widen.
///
/// ### Returns
///
/// The value as `f64`.
fn to_f64<F: BixverseFloat>(x: F) -> f64 {
    x.to_f64().unwrap_or(f64::NAN)
}

/// Widen a slice of floats to `f64`.
///
/// ### Params
///
/// * `x` - The values to widen.
///
/// ### Returns
///
/// The values as `f64`.
fn floats_to_f64<F: BixverseFloat>(x: &[F]) -> Vec<f64> {
    x.iter().map(|&v| to_f64(v)).collect()
}
