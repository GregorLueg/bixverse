use extendr_api::prelude::*;

use rustc_hash::FxHashMap;

use crate::core::enrichment::gsva::*;
use crate::utils::r_rust_interface::{faer_to_r_matrix, r_matrix_to_faer};

/// Prepare a pathway list for GSVA
///
/// @param feature_names String vector. The feature names of the matrix (should
/// be the rows).
/// @param pathway_list List. A list containing the pathways and the respective
/// genes
/// @param min_size,max_size Integer. The minimum and maximum size respectively.
///
/// @return Returns a list with (zero-indexed) indices.
///
/// @export
#[extendr]
fn rs_prepare_gsva_gs(
    feature_names: Vec<String>,
    pathway_list: List,
    min_size: usize,
    max_size: usize,
) -> extendr_api::Result<List> {
    // HashMap for fast look ups
    let gene_map: FxHashMap<&str, usize> = feature_names
        .iter()
        .enumerate()
        .map(|(i, gene)| (gene.as_str(), i))
        .collect();

    let list_names: Vec<String> = pathway_list
        .names()
        .unwrap()
        .map(|s| s.to_string())
        .collect();

    let mut filtered_pathways = Vec::new();
    let mut filtered_names = Vec::new();

    let iterator = 0..pathway_list.len();

    for i in iterator {
        let element = pathway_list.elt(i)?;
        if let Some(internal_vals) = element.as_string_vector() {
            let mut indices = Vec::with_capacity(internal_vals.len());

            for gene in &internal_vals {
                if let Some(&idx) = gene_map.get(gene.as_str()) {
                    indices.push(idx as i32);
                }
            }

            if (indices.len() >= min_size) & (indices.len() <= max_size) {
                indices.sort_unstable();
                filtered_pathways.push(indices);
                filtered_names.push(list_names[i].clone());
            }
        }
    }

    let mut result_list = List::new(filtered_pathways.len());
    for (idx, pathway) in filtered_pathways.into_iter().enumerate() {
        result_list.set_elt(idx, Robj::from(pathway))?;
    }
    result_list.set_names(filtered_names)?;

    Ok(result_list)
}

/// Rust version of the GSVA algorithm
///
/// @description
/// Rust-based implementation of the popular GSVA algorithm. Has further
/// performance optimisations compared to the original implementation.
///
/// @param exp Numerical matrix. The expression matrix with rows = genes, and
/// columns = samples
/// @param gs_list List. A list containing the indices of the pathway genes
/// (needs to be null indexed). See [bixverse::rs_prepare_gsva_gs()].
/// @param tau Float. Tau parameter, usual recommendation is to use `1.0` here.
/// Larger values emphasise the tails more.
/// @param gaussian Boolean. If `TRUE` the Gaussian kernel will be used, if
/// `FALSE` the Poisson kernel will be used.
/// @param max_diff Boolean. Scoring mode: `TRUE` = difference, `FALSE` = larger
/// absolute value
/// @param abs_rank Booelan. If `TRUE` = pos-neg, `FALSE` = pos+neg
/// @param timings Boolean. Prints timings from the algorithm.
///
/// @return Returns a matrix of gene set ES scores x samples.
///
/// @export
#[extendr]
fn rs_gsva(
    exp: RMatrix<f64>,
    gs_list: List,
    tau: f64,
    gaussian: bool,
    max_diff: bool,
    abs_rank: bool,
    timings: bool,
) -> extendr_api::Result<extendr_api::RArray<f64, [usize; 2]>> {
    let exp = r_matrix_to_faer(&exp);

    let gs_indices = get_gsva_gs_indices(gs_list)?;

    let results = gsva(
        &exp,
        &gs_indices,
        gaussian,
        tau,
        max_diff,
        abs_rank,
        timings,
    );

    Ok(faer_to_r_matrix(results.as_ref()))
}

/// Rust version of the ssGSEA algorithm
///
/// @description
/// Rust-based implementation of the popular single sample GSEA algorithm. Has
/// further performance optimisations compared to the original implementation.
///
/// @param exp Numerical matrix. The expression matrix with rows = genes, and
/// columns = samples
/// @param gs_list List. A list containing the indices of the pathway genes
/// (needs to be null indexed). See [bixverse::rs_prepare_gsva_gs()].
/// @param alpha Float. The alpha parameter to adjust the weights.
/// @param normalise Boolean. Shall the scores be normalised.
/// @param timings Boolean. Prints timings from the algorithm.
///
/// @return Returns a matrix of gene set ES scores x samples.
///
/// @export
#[extendr]
fn rs_ssgsea(
    exp: RMatrix<f64>,
    gs_list: List,
    alpha: f64,
    normalise: bool,
    timings: bool,
) -> extendr_api::Result<extendr_api::RArray<f64, [usize; 2]>> {
    let exp = r_matrix_to_faer(&exp);

    let gs_indices = get_gsva_gs_indices(gs_list)?;

    let results = ssgsea(&exp, &gs_indices, alpha, normalise, timings);

    Ok(faer_to_r_matrix(results.as_ref()))
}

extendr_module! {
    mod r_gsva;
    fn rs_prepare_gsva_gs;
    fn rs_gsva;
    fn rs_ssgsea;
}
