//! Rust to R interface for the trajectory inference methods, i.e. Palantir and
//! PAGA.
//!
//! Both algorithms take a precomputed kNN graph and nothing else, so the
//! bindings are handed the `SingleCellNearestNeighbour` list straight from R.
//! Every cell index that crosses this boundary, in and out, is 0-indexed.

use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::sc_trajectory::gene_trends::{
    compute_gene_trends, select_branch_cells, BranchSelectionParams, GeneTrendsParams,
};
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
    fn rs_gene_trends;
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
/// (0-indexed!), `dist`, `k` and `dist_metric`.
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
///  scaled to `[0, 1]`. The start cell is not pinned to 0; a start cell far
///  from 0 means the refinement disagreed with the anchor.
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
///  \item converged - Boolean. Did the pseudotime refinement converge before
///  the cap.
///  \item eigen_converged - Boolean. Did the diffusion eigensolve meet its
///  tolerance rather than running out of restarts. `FALSE` means the embedding
///  is under-resolved and every distance taken on it is suspect.
///  \item eigen_residual - Numeric. Largest achieved
///  `||A x - lambda x||` from the diffusion eigensolve.
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
    let (knn_indices, knn_distances, _, _) = knn_data_to_rust(knn_data)?;
    let params = PalantirParams::from_r_list(palantir_params)?;

    let terminal_states: Option<Vec<usize>> = match terminal_states {
        Nullable::Null => None,
        Nullable::NotNull(x) => Some(x.as_slice().r_int_convert()),
    };

    let res = run_palantir(
        &knn_indices,
        &knn_distances,
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
        eigen_converged = res.eigen_converged,
        eigen_residual = res.eigen_residual,
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

    let connectivities = sparse_data_to_list(res.connectivities).to_extendr()?;
    let connectivities_tree = sparse_data_to_list(res.connectivities_tree).to_extendr()?;

    Ok(list!(
        connectivities = connectivities,
        connectivities_tree = connectivities_tree,
        sizes = res.sizes.as_slice().r_int_convert()
    ))
}

/////////////////
// Gene trends //
/////////////////

/// Fit Palantir gene trends over pseudotime
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Selects the cells belonging to each branch from the fate probabilities, then
/// fits a landmark Gaussian process with a Matern-5/2 kernel per branch. This
/// is the Mellon-based estimator of the reference, not the legacy GAM.
///
/// Whatever expression matrix is handed over is what gets fitted. Nothing here
/// imputes; that decision belongs to the caller.
///
/// The defaults are prior-dominated. Palantir's pseudotime is min-max scaled to
/// `[0, 1]`, so a `length_scale` of `1.0` spans the whole domain and a `sigma`
/// of `1.0` sits at roughly the signal scale of log-normalised expression. That
/// resolves almost any gene into a smooth monotone or single-peaked curve.
/// Shorten `length_scale` before believing a bump.
///
/// @param expression Numerical matrix of cells x genes. All cells, in the same
/// row order as `pseudotime`, not just the branch members.
/// @param pseudotime Numerical vector. Pseudotime per cell.
/// @param branch_probs Numerical matrix of cells x fates with the fate
/// probabilities. Rows need not sum to one.
/// @param branch_params List. Parameter list, see
/// [bixverse::params_sc_branch_selection()].
/// @param gene_trend_params List. Parameter list, see
/// [bixverse::params_sc_gene_trends()].
///
/// @returns A list with
/// \itemize{
///  \item trends - List of numerical matrices, one per branch, each of
///  resolution x genes.
///  \item grids - List of numerical vectors with the pseudotime grid per
///  branch, running from the branch minimum to its maximum.
///  \item branch_cells - List of integer vectors with the cell indices
///  (0-indexed!) selected for each branch.
///  \item n_cells - Integer vector with the cell count per branch.
///  \item jitter_used - Numerical vector with the jitter each branch's
///  Cholesky needed.
/// }
///
/// @references Setty, et al., Nat. Biotechnol., 2019.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_gene_trends(
    expression: RMatrix<f64>,
    pseudotime: &[f64],
    branch_probs: RMatrix<f64>,
    branch_params: List,
    gene_trend_params: List,
) -> Result<List> {
    let selection = BranchSelectionParams::from_r_list(branch_params)?;
    let params = GeneTrendsParams::from_r_list(gene_trend_params)?;

    // the Rust side works in f32 for both, so this is the one copy we pay for
    let pseudotime: Vec<f32> = pseudotime.iter().map(|x| *x as f32).collect();
    let branch_probs = r_matrix_to_faer_fp32(&branch_probs);
    let expression = r_matrix_to_faer(&expression);

    let branch_cells =
        select_branch_cells(&pseudotime, branch_probs.as_ref(), &selection).to_extendr()?;

    let res = compute_gene_trends(
        expression,
        &pseudotime,
        &branch_cells,
        Some(branch_probs.as_ref()),
        &params,
    )
    .to_extendr()?;

    let n_branches = res.trends.len();

    let mut trends = List::new(n_branches);
    let mut grids = List::new(n_branches);
    let mut cells = List::new(n_branches);

    for i in 0..n_branches {
        trends.set_elt(i, Robj::from(faer_to_r_matrix(res.trends[i].as_ref())))?;
        grids.set_elt(i, Robj::from(res.grids[i].clone()))?;
        cells.set_elt(i, Robj::from(branch_cells[i].as_slice().r_int_convert()))?;
    }

    Ok(list!(
        trends = trends,
        grids = grids,
        branch_cells = cells,
        n_cells = res.n_cells.as_slice().r_int_convert(),
        jitter_used = res.jitter_used
    ))
}
