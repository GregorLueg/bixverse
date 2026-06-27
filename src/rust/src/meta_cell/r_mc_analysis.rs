//! Analysis methods for meta cells specifically. At the moment, it supports
//! mostly SCENIC.

use bixverse_rs::methods::nmf_hals::HalsOpts;
use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::mc_analysis::aucell::calculate_aucell_metacells;
use bixverse_rs::single_cell::mc_analysis::nmf_mc::*;
use bixverse_rs::single_cell::mc_analysis::scenic_metacells::run_scenic_grn_in_memory;
use bixverse_rs::single_cell::sc_analysis::scenic::ScenicParams;
use extendr_api::*;
use faer::Mat;

use crate::meta_cell::utils::*;

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_mc_analysis;
    // scenic
    fn rs_mc_scenic;
    // aucell
    fn rs_mc_aucell;
    // nmf
    fn rs_nmf_single_mc;
    fn rs_nmf_multi_mc;
}

/////////////////////
// Metacell SCENIC //
/////////////////////

/// SCENIC on MetaCells
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Assumes that the sparse data is pre-filtered for the genes you wish to
/// include. The indices need to be 0-indexed.
///
/// @param sparse_data A named list that needs to have `data`, `indptr`,
/// `indices`, `nrow`, `ncol` and `format`.
/// @param tf_indices Integer vector. The indices of the transcription factors.
/// @param scenic_params Named list. Contains all of the parameters need for
/// SCENIC.
/// @param seed Integer. Controls reproducibility of the function.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
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
    verbose: usize,
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
/// `r lifecycle::badge("experimental")`
/// The function will take in a list of gene set indices (0-indexed!) and
/// calculate an AUCell type statistic. Two options here: calculate this
/// with proper AUROC calculations (useful for marker gene expression) or
/// based on the Mann-Whitney statistic (useful for pathway activity
/// measurs). This version works on MetaCell counts which are stored in memory
/// directly.
///
/// @param sparse_data A named list that needs to have `data`, `indptr`,
/// `indices`, `nrow`, `ncol` and `format`.
/// @param gs_list List. List with the gene set indices (0-indexed!) of the
/// genes of interest.
/// @param auc_type String. One of `"wilcox"` or `"auroc"`, pending on
/// which statistic you wish to calculate.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
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
    verbose: usize,
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

    let res = calculate_aucell_metacells(&sparse, &gs_indices, auc_type, verbose).to_extendr()?;

    let auc_mat = Mat::from_fn(res[0].len(), res.len(), |i, j| res[j][i] as f64);
    Ok(faer_to_r_matrix(auc_mat.as_ref()))
}

/////////
// NMF //
/////////

/// Run NMF (HALS) on MetaCells
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Assumes that the sparse data is pre-filtered for the cells/genes you wish
/// to include. Indices in the sparse data need to be 0-indexed.
///
/// @param sparse_data A named list with `data`, `indptr`, `indices`, `nrow`,
/// `ncol` and `format`.
/// @param k Integer. Number of latent factors to return.
/// @param preprocessing String. One of `c("none", "sd", "sqrt_sd")`.
/// @param use_second_layer Boolean. If `TRUE`, runs NMF on normalised counts.
/// @param nmf_hals_params Named list. Contains the NMF parameters.
/// @param seed Integer. Random seed for initialisation.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with `w`, `h`, `final_loss`, `n_iter`, `converged`.
///
/// @export
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_nmf_single_mc(
    sparse_data: List,
    k: usize,
    preprocessing: &str,
    use_second_layer: bool,
    nmf_hals_params: List,
    seed: usize,
    verbose: usize,
) -> Result<List> {
    let sparse: CompressedSparseData2<f64, f64> =
        list_to_sparse_matrix(sparse_data, true).to_extendr()?;
    let sparse = cast_compressed_sparse_data_f32(sparse);
    let nmf_hals_opt: HalsOpts<f32> = HalsOpts::from_r_list(nmf_hals_params, seed).to_extendr()?;
    let nmf_res = nmf_single_run_mc(
        sparse,
        k,
        preprocessing,
        use_second_layer,
        Some(nmf_hals_opt),
        verbose,
    )
    .to_extendr()?;
    Ok(list!(
        w = faer_to_r_matrix(nmf_res.w.as_ref()),
        h = faer_to_r_matrix(nmf_res.h.as_ref()),
        final_loss = nmf_res.final_loss as f64,
        n_iter = nmf_res.n_iter,
        converged = nmf_res.converged
    ))
}

/// Run multiple NMF (HALS) restarts on MetaCells
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Assumes that the sparse data is pre-filtered for the cells/genes you wish
/// to include. Indices in the sparse data need to be 0-indexed.
///
/// @param sparse_data A named list with `data`, `indptr`, `indices`, `nrow`,
/// `ncol` and `format`.
/// @param k Integer. Number of latent factors per run.
/// @param preprocessing String. One of `c("none", "sd", "sqrt_sd")`.
/// @param use_second_layer Boolean. If `TRUE`, runs NMF on normalised counts.
/// @param nmf_hals_params Named list. Contains the NMF parameters.
/// @param n_runs Integer. Number of random restarts.
/// @param seed Integer. Base random seed. Run `i` uses `seed + i`.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with `w_all`, `h_per_run`, `losses`, `converged`,
/// `best_idx` (1-indexed).
///
/// @export
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_nmf_multi_mc(
    sparse_data: List,
    k: usize,
    preprocessing: &str,
    use_second_layer: bool,
    nmf_hals_params: List,
    n_runs: usize,
    seed: usize,
    verbose: usize,
) -> Result<List> {
    // CHECK: same caveat as rs_nmf_single_mc on the f32 conversion.
    let sparse: CompressedSparseData2<f64, f64> =
        list_to_sparse_matrix(sparse_data, true).to_extendr()?;
    let sparse = cast_compressed_sparse_data_f32(sparse);
    let nmf_hals_opt: HalsOpts<f32> = HalsOpts::from_r_list(nmf_hals_params, seed).to_extendr()?;
    let nmf_res = nmf_multiple_run_mc(
        sparse,
        k,
        preprocessing,
        use_second_layer,
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
        losses = nmf_res.losses.r_float_convert(),
        converged = nmf_res.converged,
        best_idx = (nmf_res.best_idx + 1) as i32
    ))
}
