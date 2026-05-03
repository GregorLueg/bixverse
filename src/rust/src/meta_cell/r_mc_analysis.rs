//! Analysis methods for meta cells specifically. At the moment, it supports
//! mostly SCENIC.

use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::mc_analysis::scenic_metacells::run_scenic_grn_in_memory;
use bixverse_rs::single_cell::sc_analysis::scenic::ScenicParams;
use extendr_api::*;

use crate::meta_cell::utils::cast_compressed_sparse_data_u32;

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_mc_analysis;
    // scenic
    fn rs_mc_scenic;
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
/// @param gene_indices Integer vector. If you want to run the algorithm on a
/// subset of cells.
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
    gene_indices: Vec<i32>,
    tf_indices: Vec<i32>,
    scenic_params: List,
    seed: usize,
    verbose: bool,
) -> Result<RArray<f64, 2>> {
    let gene_indices = gene_indices.r_int_convert();
    let tf_indices = tf_indices.r_int_convert();
    let sparse: CompressedSparseData2<f64, f64> =
        list_to_sparse_matrix(sparse_data, true).to_extendr()?;
    let sparse = cast_compressed_sparse_data_u32(sparse);
    let scenic_params = ScenicParams::from_r_list(scenic_params)?;

    let grn_matrix = run_scenic_grn_in_memory(
        &sparse,
        &gene_indices,
        &tf_indices,
        &scenic_params,
        seed,
        verbose,
    )
    .to_extendr()?;

    Ok(faer_to_r_matrix(grn_matrix.as_ref()))
}
