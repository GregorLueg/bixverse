//! Functions related to cell type label transfer or cell type annotations

use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::sc_annotation::sc_type::*;
use extendr_api::*;

use crate::single_cell::utils::process_cell_markers;

////////////////////
// extendr Module //
////////////////////

extendr_module! {
    // module
    mod r_sc_annotation;
    // sctype
    fn rs_sc_type;
    fn rs_sc_type_cluster_assignment;
}

////////////
// ScType //
////////////

/// Run the ScType scoring approach
///
/// @description This Rust function implements the cell type scoring approach
/// from Ianevski et al. (2022).
///
/// @param f_path String. Path to the `counts_genes.bin` file.
/// @param cell_indices Integer vector. 0-indexed(!) positions of cells to
/// include in the analysis
/// @param cell_markers A list with the cell marker gene indices.
/// @param sensitivity Boolean. Shall a sensitivity correction be applied that
/// downweights common cell type markers.
/// @param weight_floor Optional numeric. If `sensitivity = TRUE`, what is
/// the weight floor. If not provided, defaults to `0.1`.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with
/// \itemize{
///   \item cell_types - String vector. The cell types
///   \item scores - Row-major scores (cells x cell_types).
///   \item n_cells - Number of cells
///   \item n_cell_types - Number of cell types
/// }
///
/// @export
#[extendr]
fn rs_sc_type(
    f_path: &str,
    cell_indices: Vec<i32>,
    cell_markers: List,
    sensitivity: bool,
    weight_floor: Option<f64>,
    verbose: usize,
) -> Result<List> {
    let cell_markers = process_cell_markers(cell_markers)?;
    let cell_indices = cell_indices.r_int_convert();

    let res: SctypeRes = run_sctype(
        f_path,
        &cell_indices,
        &cell_markers,
        sensitivity,
        None,
        weight_floor.map(|x| x as f32),
        verbose,
    )
    .to_extendr()?;

    Ok(list!(
        cell_types = res.cell_types,
        scores = res.scores.r_float_convert(),
        n_cells = res.n_cells,
        n_cell_types = res.n_cell_types
    ))
}

/// Score the individual clusters based on ScType
///
/// @description This Rust function implements the cell type scoring approach
/// from Ianevski et al. (2022).
///
/// @param sc_type_res List. The ScType results.
/// @param cluster_labels Integer. Cluster assignment. Needs to be of length
/// of scored cells.
///
/// @returns A list with
/// \itemize{
///  \item cluster_id - The cluster id/integer
///  \item cell_type - String; the predicted cell type
///  \item score - The final score for the clsuter.
///  \item n_cells - The number of cells in the cluster.
/// }
#[extendr]
fn rs_sc_type_cluster_assignment(sc_type_res: List, cluster_labels: Vec<i32>) -> Result<List> {
    let cluster_labels = cluster_labels.r_int_convert();

    let sc_res = SctypeRes::from_r_list(sc_type_res)?;

    let cluster_assignments: Vec<ScTypeClusterAssignment> =
        assign_clusters(&sc_res, &cluster_labels).to_extendr()?;

    let mut cluster_id: Vec<usize> = Vec::with_capacity(cluster_assignments.len());
    let mut cell_type: Vec<String> = Vec::with_capacity(cluster_assignments.len());
    let mut scores: Vec<f64> = Vec::with_capacity(cluster_assignments.len());
    let mut n_cells: Vec<usize> = Vec::with_capacity(cluster_assignments.len());

    for res in cluster_assignments {
        cluster_id.push(res.cluster);
        cell_type.push(res.cell_type);
        scores.push(res.score as f64);
        n_cells.push(res.n_cells);
    }

    Ok(list!(
        cluster_id = cluster_id.r_int_convert(),
        cell_type = cell_type,
        scores = scores,
        n_cells = n_cells.r_int_convert()
    ))
}
