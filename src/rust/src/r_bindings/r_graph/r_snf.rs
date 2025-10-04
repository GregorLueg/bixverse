use extendr_api::scalar::Rbool;
use extendr_api::*;
use faer::MatRef;

use crate::core::graph::snf::*;
use crate::core::methods::rbh::r_matrix_list_to_vec;
use crate::utils::r_rust_interface::*;

/// Calculate the SNF affinity matrix for continuous values
///
/// @param data Numerical matrix. Needs to be oriented features x samples!
/// @param distance_type String. One of
/// `c("euclidean", "manhattan", "canberra", "cosine")`. Which distance metric
/// to use here.
/// @param k Integer. Number of neighbours to consider.
/// @param mu Float. Normalisation factor for the Gaussian kernel width.
/// @param normalise Boolean. Shall continuous values be Z-scored.
///
/// @return The affinity matrix based on continuous values.
///
/// @export
#[extendr]
fn rs_snf_affinity_continuous(
    data: RMatrix<f64>,
    distance_type: &str,
    k: usize,
    mu: f64,
    normalise: bool,
) -> extendr_api::Result<RArray<f64, [usize; 2]>> {
    let data = r_matrix_to_faer(&data);

    let affinity_data = make_affinity_continuous(&data, distance_type, k, mu, normalise)?;

    Ok(faer_to_r_matrix(affinity_data.as_ref()))
}

/// Calculate the SNF affinity matrix for categorical values
///
/// @param data Integer matrix. Needs to be oriented features x samples! The
/// integers represent the factor values of the catagories.
/// @param k Integer. Number of neighbours to consider.
/// @param mu Float. Normalisation factor for the Gaussian kernel width.
///
/// @return The affinity matrix based on categorical values.
///
/// @export
#[extendr]
fn rs_snf_affinity_cat(
    data: RMatrix<i32>,
    k: usize,
    mu: f64,
) -> extendr_api::Result<RArray<f64, [usize; 2]>> {
    let data = r_matrix_to_faer(&data);

    let affinity_data = make_affinity_categorical(&data, k, mu)?;

    Ok(faer_to_r_matrix(affinity_data.as_ref()))
}

/// Calculate the SNF affinity matrix for mixed values
///
/// @param data Numerical matrix. Needs to be oriented features x samples! This
/// function will calculate the Gower distance under the hood for the affinity
/// calculation.
/// @param is_cat Boolean vector. Which of the features are categorical. Needs
/// to be of `nrow(data)`.
/// @param k Integer. Number of neighbours to consider.
/// @param mu Float. Normalisation factor for the Gaussian kernel width.
///
/// @return The affinity matrix based on mixed values.
///
/// @export
#[extendr]
fn rs_snf_affinity_mixed(
    data: RMatrix<f64>,
    is_cat: Vec<Rbool>,
    k: usize,
    mu: f64,
) -> extendr_api::Result<RArray<f64, [usize; 2]>> {
    let data = r_matrix_to_faer(&data);

    let is_cat = is_cat
        .iter()
        .map(|r_obj| r_obj.to_bool())
        .collect::<Vec<bool>>();

    let affinity_data = make_affinity_mixed(&data, &is_cat, k, mu)?;

    Ok(faer_to_r_matrix(affinity_data.as_ref()))
}

/// Similarity network fusion
///
/// @description This function iteratively fuses the affinity matrices together.
///
/// @param aff_mat_list A list of numerical matrices. The affinity matrices to
/// fuse together.
/// @param k Integer. Number of neighbours to consider.
/// @param t Integer. Number of iterations for the algorithm.
/// @param alpha Float. Normalisation parameter controlling the fusion strength.
///
/// @return The final affinity matrix after the fusion.
///
/// @export
#[extendr]
fn rs_snf(
    aff_mat_list: List,
    k: usize,
    t: usize,
    alpha: f64,
) -> extendr_api::Result<RArray<f64, [usize; 2]>> {
    let mat_list = r_matrix_list_to_vec(aff_mat_list);

    // Store owned matrices to keep them alive
    let owned_mats: Vec<_> = mat_list.into_iter().map(|(_, mat)| mat).collect();

    // Create references from the owned matrices
    let aff_mats: Vec<MatRef<f64>> = owned_mats.iter().map(r_matrix_to_faer).collect();

    let res = snf(&aff_mats, k, t, alpha)?;

    Ok(faer_to_r_matrix(res.as_ref()))
}

extendr_module! {
    mod r_snf;
    fn rs_snf_affinity_continuous;
    fn rs_snf_affinity_cat;
    fn rs_snf_affinity_mixed;
    fn rs_snf;
}
