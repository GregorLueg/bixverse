//! Analysis methods for meta cells specifically. At the moment, it supports
//! mostly SCENIC.

use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::mc_analysis::aucell::calculate_aucell_metacells;
use bixverse_rs::single_cell::mc_analysis::scenic_metacells::run_scenic_grn_in_memory;
use bixverse_rs::single_cell::sc_analysis::scenic::ScenicParams;
use extendr_api::*;
use faer::Mat;

use crate::meta_cell::utils::cast_compressed_sparse_data_u32;

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_mc_analysis;
    // scenic
    fn rs_mc_scenic;
    // aucell
    fn rs_mc_aucell;
}

/////////////////////
// Metacell SCENIC //
/////////////////////

/// SCENIC on MetaCells
///
/// @description
/// Assumes that the sparse data is pre-filtered for the genes you wish to
/// include. The indices need to be 0-indexed.
///
/// @param sparse_data A named list that needs to have `data`, `indptr`,
/// `indices`, `nrow`, `ncol` and `format`.
/// @param tf_indices Integer vector. The indices of the transcription factors.
/// @param scenic_params Named list. Contains all of the parameters need for
/// SCENIC.
/// @param seed Integer. Controls reproducibility of the function.
/// @param verbose Boolean. Controls the verbosity of the function.
///
/// @returns A gene x TF importance matrix
///
/// @export
#[extendr]
fn rs_mc_scenic(
    sparse_data: List,
    tf_indices: Vec<i32>,
    scenic_params: List,
    seed: usize,
    verbose: bool,
) -> Result<RArray<f64, 2>> {
    let tf_indices = tf_indices.r_int_convert();
    let sparse: CompressedSparseData2<f64, f64> =
        list_to_sparse_matrix(sparse_data, true).to_extendr()?;
    let sparse = cast_compressed_sparse_data_u32(sparse);
    let scenic_params = ScenicParams::from_r_list(scenic_params)?;

    let grn_matrix = run_scenic_grn_in_memory(&sparse, &tf_indices, &scenic_params, seed, verbose)
        .to_extendr()?;

    Ok(faer_to_r_matrix(grn_matrix.as_ref()))
}

/////////////////////
// Metacell AUCell //
/////////////////////

/// Calculate AUCell in Rust (for meta cells)
///
/// @description
/// The function will take in a list of gene set indices (0-indexed!) and
/// calculate an AUCell type statistic. Two options here: calculate this
/// with proper AUROC calculations (useful for marker gene expression) or
/// based on the Mann-Whitney statistic (useful for pathway activity
/// measurs). This version works on MetaCell counts which are stored in memory
/// directly.
///
/// @param gs_list List. List with the gene set indices (0-indexed!) of the
/// genes of interest.
/// @param cells_to_keep Integer. Vector of indices of the cells to keep.
/// @param auc_type String. One of `"wilcox"` or `"auroc"`, pending on
/// which statistic you wish to calculate.
/// @param streaming Boolean. Shall the data be streamed.
/// @param verbose Boolean. Controls verbosity of the function.
///
/// @return A matrix of cells x gene sets with the values representing the
/// AUC.
///
/// @export
#[extendr]
fn rs_mc_aucell(
    sparse_data: List,
    gs_list: List,
    auc_type: &str,
    verbose: bool,
) -> Result<RArray<f64, 2>> {
    let mut gs_indices: Vec<Vec<usize>> = Vec::with_capacity(gs_list.len());
    for i in 0..gs_list.len() {
        let r_obj = gs_list.elt(i).unwrap();
        let int = r_obj
            .as_integer_vector()
            .unwrap()
            .iter()
            .map(|x| *x as usize)
            .collect::<Vec<usize>>();
        gs_indices.push(int);
    }

    let sparse: CompressedSparseData2<f64, f64> =
        list_to_sparse_matrix(sparse_data, true).to_extendr()?;
    let sparse = cast_compressed_sparse_data_u32(sparse);

    let res = calculate_aucell_metacells(&sparse, &gs_indices, auc_type, verbose);

    let auc_mat = Mat::from_fn(res[0].len(), res.len(), |i, j| res[j][i] as f64);
    Ok(faer_to_r_matrix(auc_mat.as_ref()))
}
