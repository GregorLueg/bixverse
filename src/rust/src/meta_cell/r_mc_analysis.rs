//! Analysis methods for meta cells specifically. At the moment, it supports
//! mostly SCENIC.

use bixverse_rs::methods::nmf_hals::HalsOpts;
use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::mc_analysis::aucell::calculate_aucell_metacells;
use bixverse_rs::single_cell::mc_analysis::hotspot_mc::*;
use bixverse_rs::single_cell::mc_analysis::nmf_mc::*;
use bixverse_rs::single_cell::mc_analysis::scenic_metacells::run_scenic_grn_in_memory;
use bixverse_rs::single_cell::sc_analysis::dge_pathway_scores::AucellParams;
use bixverse_rs::single_cell::sc_analysis::hotspot::{HotSpotGeneRes, HotSpotPairRes, HotSpotParams};
use bixverse_rs::single_cell::sc_analysis::scenic::ScenicParams;
use extendr_api::*;
use faer::{Mat, MatRef};

use crate::meta_cell::utils::*;
use crate::single_cell::utils::knn_data_to_rust;

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_mc_analysis;
    // scenic
    fn rs_mc_scenic;
    // aucell
    fn rs_mc_aucell;
    // hotspot
    fn rs_mc_hotspot_autocor;
    fn rs_mc_hotspot_gene_cor;
    // nmf
    fn rs_nmf_single_mc;
    fn rs_nmf_multi_mc;
}

/////////////
// Helpers //
/////////////

/// Type for a resolved kNN graph
///
/// ### Fields
///
/// * `0` - Neighbour indices per meta cell, self excluded
/// * `1` - The matching neighbour distances, ascending
/// * `2` - Whether those distances hold `d^2`
type ResolvedKnn = extendr_api::Result<(Vec<Vec<usize>>, Vec<Vec<f32>>, bool)>;

/// Resolve the kNN graph the HotSpot kernel runs over.
///
/// Either the pre-computed graph handed over from R, or one built from the
/// embedding. Whether the distances are pre-squared follows from the metric,
/// which is read off the supplied graph rather than the parameter list: a
/// cached graph may well have been built with a different metric than
/// `ann_dist` says.
///
/// ### Params
///
/// * `knn_data` - Optional pre-computed kNN data with `indices`, `dist` and
///   `dist_metric`
/// * `embd` - The embedding, metacells x dimensions
/// * `knn_params` - Parameters for the ANN search, only read when a graph is
///   built
/// * `seed` - Random seed for the ANN search
/// * `verbosity` - Parsed verbosity level
///
/// ### Returns
///
/// The neighbour indices, their distances, and whether those distances hold
/// `d^2`.
fn resolve_knn_graph(
    knn_data: Nullable<List>,
    embd: MatRef<f32>,
    knn_params: &KnnParams,
    seed: usize,
    verbosity: Verbosity,
) -> ResolvedKnn {
    if knn_data != extendr_api::Nullable::Null {
        if verbosity.normal_verbosity() {
            println!("Using provided kNN graph...")
        }
        let knn_data = knn_data
            .into_robj()
            .as_list()
            .ok_or_else(|| Error::Other("'knn_data' is not a list".into()))?;
        let (knn_indices, knn_dist, _, dist_metric) = knn_data_to_rust(knn_data)?;

        Ok((knn_indices, knn_dist, distances_are_squared(&dist_metric)))
    } else {
        if verbosity.normal_verbosity() {
            println!("Generating a kNN graph from scratch")
        }
        let (knn_indices, knn_dist) = generate_knn_with_dist(
            embd,
            knn_params,
            true,
            false,
            seed,
            verbosity.detailed_verbosity(),
        )
        .to_extendr()?;

        Ok((
            knn_indices,
            knn_dist.unwrap(),
            distances_are_squared(&knn_params.ann_dist),
        ))
    }
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
/// calculate an AUCell type statistic. Three options here: the recovery-curve
/// AUC of Aibar, et al. (the actual AUCell statistic), an AUC derived from the
/// Mann-Whitney statistic, or average precision. This version works on
/// MetaCell counts which are stored in memory directly.
///
/// @param sparse_data A named list that needs to have `data`, `indptr`,
/// `indices`, `nrow`, `ncol` and `format`.
/// @param gs_list List. List with the gene set indices (0-indexed!) of the
/// genes of interest.
/// @param aucell_params List. The AUCell parameters, see
/// [bixverse::params_sc_aucell()].
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
    aucell_params: List,
    verbose: usize,
) -> Result<RArray<f64, 2>> {
    let aucell_params = AucellParams::from_r_list(aucell_params)?;

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

    let res = calculate_aucell_metacells(&sparse, &gs_indices, Some(aucell_params), verbose)
        .to_extendr()?;

    let auc_mat = Mat::from_fn(res[0].len(), res.len(), |i, j| res[j][i] as f64);
    Ok(faer_to_r_matrix(auc_mat.as_ref()))
}

//////////////////////
// Metacell HotSpot //
//////////////////////

/// Calculate gene spatial auto-correlations (for meta cells)
///
/// @description
/// `r lifecycle::badge("experimental")`
/// This function implements the HotSpot auto-correlation functionality and
/// will return to what extent a given gene shows auto-correlation in the
/// kNN-graph over the meta cells. For details see DeTomaso, et al. This version
/// works on MetaCell counts which are stored in memory directly. There is no
/// streaming variant: streaming bounds disk re-reads, which is not a problem
/// an in-memory matrix has.
///
/// @param sparse_data A named list that needs to have `data`, `indptr`,
/// `indices`, `nrow`, `ncol` and `format`. Shape is (metacells, genes) and the
/// data are the raw counts.
/// @param embd Numerical matrix. The embedding matrix from which to generate
/// the kNN graph.
/// @param knn_data Optional list. This contains pre-computed kNN data
/// (including distances) and the `dist_metric` it was built with. The user has
/// to ensure consistency! If provided, this will be used and whether the
/// distances are treated as squared is derived from `dist_metric` rather than
/// from the parameter list.
/// @param hotspot_params List. The HotSpot parameter list. The kNN parameters
/// are only read when no `knn_data` is provided.
/// @param cells_to_keep Integer vector. 0-index vector indicating which meta
/// cells to include in the analysis. Ensure that this is of same order/length
/// as the embedding matrix.
/// @param genes_to_use Integer vector. 0-index vector indicating which genes
/// to include.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
/// @param seed Integer. Random seed for reproducibility.
///
/// @returns A list with the following elements.
/// \itemize{
///   \item gene_idx - 0-based integer indicating the gene index.
///   \item gaerys_c - Gaery's C calculation for the autocorrelation
///   coefficient.
///   \item z_score - Z-score of the auto-correlation.
///   \item pval - P-value derived from the Z-score.
///   \item fdr - False discovery rate based on the p-value.
/// }
///
/// @export
///
/// @references DeTomaso, et al., Cell Systems, 2021
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_mc_hotspot_autocor(
    sparse_data: List,
    embd: RMatrix<f64>,
    knn_data: Nullable<List>,
    hotspot_params: List,
    cells_to_keep: Vec<i32>,
    genes_to_use: Vec<i32>,
    verbose: usize,
    seed: usize,
) -> extendr_api::Result<List> {
    assert!(
        embd.nrows() == cells_to_keep.len(),
        "The embedding matrix need to have the same nrow as the the cells to use."
    );

    let verbosity = parse_verbosity_level(verbose);
    let mut hotspot_params = HotSpotParams::from_r_list(hotspot_params)?;

    let embd = r_matrix_to_faer_fp32(&embd);
    let cells_to_keep = cells_to_keep.r_int_convert();
    let genes_to_use = genes_to_use.r_int_convert();

    let (knn_indices, knn_dist, squared_distances) = resolve_knn_graph(
        knn_data,
        embd.as_ref(),
        &hotspot_params.knn_params,
        seed,
        verbosity,
    )?;

    hotspot_params.graph_params.squared_distances = squared_distances;

    // HotSpot only ever reads the raw counts and the per-cell library sizes, so
    // the second layer the reader insists on can stay the copy of the raw one
    // `list_to_sparse_matrix` puts there
    let sparse: CompressedSparseData2<f64, f64> =
        list_to_sparse_matrix(sparse_data, true).to_extendr()?;
    let sparse = cast_compressed_sparse_data_u32(sparse);

    let res: HotSpotGeneRes = hotspot_autocor_metacells(
        &sparse,
        &knn_indices,
        &knn_dist,
        &cells_to_keep,
        &genes_to_use,
        &hotspot_params,
        verbose,
    )
    .to_extendr()?;

    Ok(list!(
        gene_idx = res.gene_idx,
        gaerys_c = res.c,
        z_score = res.z,
        pval = res.pval,
        fdr = res.fdr
    ))
}

/// Calculate gene to gene spatial correlations (for meta cells)
///
/// @description
/// `r lifecycle::badge("experimental")`
/// This function implements the HotSpot gene <> gene local correlation
/// functionality from HotSpot, see DeTomaso, et al. This version works on
/// MetaCell counts which are stored in memory directly.
///
/// Three dense metacells x genes blocks are live at once, so keep
/// `genes_to_use` to the panel actually of interest rather than the whole
/// transcriptome.
///
/// @param sparse_data A named list that needs to have `data`, `indptr`,
/// `indices`, `nrow`, `ncol` and `format`. Shape is (metacells, genes) and the
/// data are the raw counts.
/// @param embd Numerical matrix. The embedding matrix from which to generate
/// the kNN graph.
/// @param knn_data Optional list. This contains pre-computed kNN data
/// (including distances) and the `dist_metric` it was built with. The user has
/// to ensure consistency! If provided, this will be used and whether the
/// distances are treated as squared is derived from `dist_metric` rather than
/// from the parameter list.
/// @param hotspot_params List. The HotSpot parameter list. The kNN parameters
/// are only read when no `knn_data` is provided; `normalise` is unused on this
/// path.
/// @param cells_to_keep Integer vector. 0-index vector indicating which meta
/// cells to include in the analysis. Ensure that this is of same order/length
/// as the embedding matrix.
/// @param genes_to_use Integer vector. 0-index vector indicating which genes
/// to include.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
/// @param seed Integer. Random seed for reproducibility.
///
/// @returns A list with the following elements.
/// \itemize{
///   \item cor - The gene x gene local correlation matrix.
///   \item z - The Z-scores of these local correlations.
/// }
///
/// @export
///
/// @references DeTomaso, et al., Cell Systems, 2021
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_mc_hotspot_gene_cor(
    sparse_data: List,
    embd: RMatrix<f64>,
    knn_data: Nullable<List>,
    hotspot_params: List,
    cells_to_keep: Vec<i32>,
    genes_to_use: Vec<i32>,
    verbose: usize,
    seed: usize,
) -> extendr_api::Result<List> {
    assert!(
        embd.nrows() == cells_to_keep.len(),
        "The embedding matrix need to have the same nrow as the the cells to use."
    );

    let verbosity = parse_verbosity_level(verbose);
    let mut hotspot_params = HotSpotParams::from_r_list(hotspot_params)?;

    let embd = r_matrix_to_faer_fp32(&embd);
    let cells_to_keep = cells_to_keep.r_int_convert();
    let genes_to_use = genes_to_use.r_int_convert();

    let (knn_indices, knn_dist, squared_distances) = resolve_knn_graph(
        knn_data,
        embd.as_ref(),
        &hotspot_params.knn_params,
        seed,
        verbosity,
    )?;

    hotspot_params.graph_params.squared_distances = squared_distances;

    let sparse: CompressedSparseData2<f64, f64> =
        list_to_sparse_matrix(sparse_data, true).to_extendr()?;
    let sparse = cast_compressed_sparse_data_u32(sparse);

    let res: HotSpotPairRes = hotspot_gene_cor_metacells(
        &sparse,
        &knn_indices,
        &knn_dist,
        &cells_to_keep,
        &genes_to_use,
        &hotspot_params,
        verbose,
    )
    .to_extendr()?;

    Ok(list!(
        cor = faer_to_r_matrix(res.cor.as_ref()),
        z = faer_to_r_matrix(res.z_scores.as_ref())
    ))
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
