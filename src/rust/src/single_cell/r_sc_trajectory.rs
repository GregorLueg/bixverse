//! Rust to R interface for the trajectory inference methods, i.e. Palantir and
//! PAGA.
//!
//! Both algorithms take a precomputed kNN graph and nothing else, so the
//! bindings are handed the `SingleCellNearestNeighbour` list straight from R.
//! Every cell index that crosses this boundary, in and out, is 0-indexed.

use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::sc_trajectory::paga::{run_paga, PagaResult};
use bixverse_rs::single_cell::sc_trajectory::palantir::{run_palantir, PalantirParams};
use extendr_api::*;

use crate::single_cell::utils::{knn_data_to_rust, knn_indices_processing};

////////////////////
// extendr Module //
////////////////////

extendr_module! {
    mod r_sc_trajectory;
    fn rs_palantir;
    fn rs_paga;
}

//////////////
// Palantir //
//////////////

/// Run Palantir trajectory inference over a kNN graph
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Implementation of Palantir in Rust. The provided kNN graph feeds the
/// diffusion kernel. The geodesics themselves are measured over a second kNN
/// graph that is built internally on the multiscale space, which is where the
/// reference measures them.
///
/// @param knn_data List. The `SingleCellNearestNeighbour` data with `indices`
/// (0-indexed!), `dist`, `k` and `dist_metric`. Whether the distances are
/// treated as squared is derived from `dist_metric`.
/// @param palantir_params List. Parameter list, see
/// [bixverse::params_sc_palantir()].
/// @param early_cell Integer. Index (0-indexed!) of the early cell within the
/// rows of the kNN data.
/// @param terminal_states Optional integer vector. Indices (0-indexed!) of the
/// terminal states. If `NULL`, they are detected from the waypoint Markov
/// chain.
/// @param seed Integer. Seed for reproducibility purposes.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with
/// \itemize{
///  \item pseudotime - Numerical vector with the pseudotime per cell, min-max
///  scaled to `[0, 1]` with the start cell at 0.
///  \item entropy - Numerical vector with the differentiation entropy per cell
///  (natural log).
///  \item branch_probs - Numerical matrix of cells x terminal states with the
///  fate probabilities. Rows need not sum to one, as sub-threshold values are
///  zeroed without renormalisation.
///  \item terminal_states - Integer vector with the terminal state cell indices
///  (0-indexed!). Sets the column order of `branch_probs`.
///  \item waypoints - Integer vector with the waypoint cell indices
///  (0-indexed!). The first element is the start cell.
///  \item start_cell - Integer. The start cell that was actually used
///  (0-indexed!).
///  \item multiscale - Numerical matrix of cells x components with the
///  multiscale diffusion components.
///  \item iterations - Integer. Refinement passes that were run.
///  \item converged - Boolean. Did the refinement converge before the cap.
///  \item repair_edges - Integer. Bridging edges the connectivity repair had to
///  add. Anything non-zero means the kNN graph was disconnected.
///  \item stranded_waypoints - Integer. Waypoints from which no terminal state
///  is reachable.
/// }
///
/// @export
///
/// @references Setty, et al., Nat. Biotechnol., 2019.
///
/// @keywords internal
#[extendr]
fn rs_palantir(
    knn_data: List,
    palantir_params: List,
    early_cell: usize,
    terminal_states: Nullable<Vec<i32>>,
    seed: usize,
    verbose: usize,
) -> Result<List> {
    let (knn_indices, knn_distances, _, distance) = knn_data_to_rust(knn_data)?;
    let params = PalantirParams::from_r_list(palantir_params)?;
    let squared_dist = distance == "euclidean";

    let terminal_states: Option<Vec<usize>> = match terminal_states {
        Nullable::Null => None,
        Nullable::NotNull(x) => Some(x.as_slice().r_int_convert()),
    };

    let res = run_palantir(
        &knn_indices,
        &knn_distances,
        squared_dist,
        early_cell,
        terminal_states.as_deref(),
        params,
        seed as u64,
        verbose,
    )
    .to_extendr()?;

    Ok(list!(
        pseudotime = res.pseudotime.r_float_convert(),
        entropy = res.entropy.r_float_convert(),
        branch_probs = faer_to_r_matrix(res.branch_probs.as_ref()),
        terminal_states = res.terminal_states.as_slice().r_int_convert(),
        waypoints = res.waypoints.as_slice().r_int_convert(),
        start_cell = res.start_cell as i32,
        multiscale = faer_to_r_matrix(res.multiscale.as_ref()),
        iterations = res.iterations as i32,
        converged = res.converged,
        repair_edges = res.repair_edges as i32,
        stranded_waypoints = res.stranded_waypoints as i32
    ))
}

//////////
// PAGA //
//////////

/// Run PAGA over a kNN graph and a clustering
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Implementation of PAGA in Rust. The kNN graph is read as a directed graph,
/// with an edge running from cell `i` to each of its neighbours. Distances are
/// not used, hence only the kNN index matrix is needed here.
///
/// @param knn_mat Integer matrix. Rows represent cells and the columns
/// represent the neighbours (0-indexed!).
/// @param partitions Integer vector. Cluster label per cell (0-indexed!). Needs
/// one entry per row of `knn_mat`.
/// @param n_partitions Integer. Declared cluster count. Pass the number of
/// factor levels to retain empty ones.
///
/// @returns A list with
/// \itemize{
///  \item connectivities - The abstracted graph, clusters x clusters, as a
///  symmetric CSR list with a zero diagonal. Use
///  [bixverse::sparse_list_to_mat()] to transform it into a sparse matrix.
///  \item connectivities_tree - Maximum spanning forest of `connectivities`,
///  carrying the original connectivity values on the retained edges. Same
///  format as above.
///  \item sizes - Integer vector with the number of cells per cluster.
/// }
///
/// @export
///
/// @references Wolf, et al., Genome Biol., 2019.
///
/// @keywords internal
#[extendr]
fn rs_paga(knn_mat: RMatrix<i32>, partitions: Vec<i32>, n_partitions: usize) -> Result<List> {
    let knn_indices = knn_indices_processing(knn_mat);
    let partitions: Vec<usize> = partitions.as_slice().r_int_convert();

    let res: PagaResult<f64> =
        run_paga(&knn_indices, &partitions, Some(n_partitions)).to_extendr()?;

    Ok(list!(
        connectivities = sparse_data_to_list(res.connectivities),
        connectivities_tree = sparse_data_to_list(res.connectivities_tree),
        sizes = res.sizes.as_slice().r_int_convert()
    ))
}
