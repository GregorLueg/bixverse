use extendr_api::prelude::*;

use faer::Mat;
use rand::prelude::*;
use rayon::prelude::*;
use std::collections::HashMap;

use crate::utils_r_rust::{faer_to_r_matrix, r_matrix_to_faer};
use crate::utils_rust::array_f64_max_min;
use crate::utils_stats::*;

/// Calculate the OT harmonic sum
///
/// @param x The numeric vector (should be between 0 and 1) for which to
/// calculate the harmonic sum
///
/// @return Returns the harmonic sum according to the OT calculation.
///
/// @export
#[extendr]
fn rs_ot_harmonic_sum(mut x: Vec<f64>) -> f64 {
    x.sort_by(|a, b| b.partial_cmp(a).unwrap());

    let harmonic_sum: f64 = x
        .iter()
        .enumerate()
        .map(|(i, x)| x / (i + 1).pow(2) as f64)
        .sum();

    let max_sum: f64 = vec![1; x.len()]
        .into_iter()
        .enumerate()
        .map(|(i, x)| x as f64 / (i + 1).pow(2) as f64)
        .sum();

    harmonic_sum / max_sum
}

/// Reconstruct a matrix from a flattened upper triangle vector
///
/// @description This function takes a flattened vector of the upper triangle
/// from a symmetric matrix (think correlation matrix) and reconstructs the full
/// dense matrix for you.
///
/// @param cor_vector Numeric vector. The vector of correlation coefficients
/// that you want to use to go back to a dense matrix.
/// @param shift Integer. If you applied a shift, i.e. included the diagonal
/// values = 0; or excluded the diagonal values = 1.
/// @param n Integer. Original dimension (i.e., ncol/nrow) of the matrix to be
/// reconstructed.
///
/// @return The dense R matrix.
///
/// @export
#[extendr]
fn rs_upper_triangle_to_dense(
    cor_vector: &[f64],
    shift: usize,
    n: usize,
) -> extendr_api::RArray<f64, [usize; 2]> {
    let mut mat = Mat::<f64>::zeros(n, n);
    let mut idx = 0;
    for i in 0..n {
        for j in i..n {
            if shift == 1 && i == j {
                mat[(i, j)] = 1_f64
            } else {
                mat[(i, j)] = cor_vector[idx];
                mat[(j, i)] = cor_vector[idx];
                idx += 1;
            }
        }
    }

    faer_to_r_matrix(mat.as_ref())
}

/// Helper function to assess CoReMo cluster quality
///
/// @description This function assesses the quality of the clusters
/// with a given cut `k`. Returns the median R2 (cor^2) and the median absolute
/// deviation (MAD) of the clusters. Large clusters (≥1000) are subsampled
/// to a random set of 1000 genes.
///
/// @param cluster_genes A list. Contains the cluster and their respective genes.
/// @param cor_mat Numerical matrix. Contains the correlation coefficients.
/// @param row_names String vector. The row names (or column names) of the
/// correlation matrix.
/// @param seed Integer. Random seed for the sub sampling of genes.
///
/// @return A list containing:
///  \itemize{
///   \item r2med - median R2 of the cluster.
///   \item r2mad - median absolute deviation of the R2 in the cluster.
///   \item size - size of the cluster.
/// }
#[extendr]
fn rs_coremo_quality(
    cluster_genes: List,
    cor_mat: RMatrix<f64>,
    row_names: Vec<String>,
    seed: u64,
) -> extendr_api::Result<List> {
    // R matrix to faer
    let cor_mat = r_matrix_to_faer(&cor_mat);

    // Faster look-ups
    let gene_map: HashMap<&str, usize> = row_names
        .iter()
        .enumerate()
        .map(|(i, gene)| (gene.as_str(), i))
        .collect();

    // Get the cluster indices
    let mut all_indices = Vec::with_capacity(cluster_genes.len());

    for i in 0..cluster_genes.len() {
        let element_i = cluster_genes.elt(i)?;
        if let Some(internal_vals) = element_i.as_string_vector() {
            let mut indices = Vec::with_capacity(internal_vals.len());
            for gene in &internal_vals {
                if let Some(&idx) = gene_map.get(gene.as_str()) {
                    indices.push(idx);
                }
            }
            all_indices.push(indices)
        }
    }

    let all_sizes: Vec<usize> = all_indices.iter().map(|vec| vec.len()).collect();

    // Iterate through the clusters and get the data
    let res: Vec<(f64, f64)> = all_indices
        .iter()
        .map(|index_vec| {
            let n = index_vec.len();

            let indices: Vec<usize> = if n > 1000 {
                let mut rng = StdRng::seed_from_u64(seed);
                index_vec.choose_multiple(&mut rng, 1000).cloned().collect()
            } else {
                index_vec.to_vec()
            };

            let indices_diagonal: Vec<(usize, &[usize])> = indices
                .iter()
                .enumerate()
                .map(|(i, first)| (*first, &index_vec[i + 1..]))
                .take_while(|(_, rest)| !rest.is_empty())
                .collect();

            let expected_size = index_vec.len() * (index_vec.len() - 1) / 2;
            let mut vals: Vec<f64> = Vec::with_capacity(expected_size);

            for (r_idx, col_idx) in indices_diagonal {
                for col in col_idx {
                    vals.push(cor_mat[(r_idx, *col)].powi(2))
                }
            }

            let r2_med = median(&vals);
            let r2_mad = mad(&vals);

            (r2_med, r2_mad)
        })
        .collect();

    let mut r2_med_vec: Vec<f64> = Vec::with_capacity(res.len());
    let mut r2_mad_vec: Vec<f64> = Vec::with_capacity(res.len());

    for (med, mad) in res {
        r2_med_vec.push(med);
        r2_mad_vec.push(mad);
    }

    Ok(list!(
        r2med = r2_med_vec,
        r2mad = r2_mad_vec,
        size = all_sizes
    ))
}

/// Apply a Radial Basis Function
///
/// @description Applies a radial basis function (RBF) to a given distance
/// vector. Has at the option to apply a Gaussian, Bump or Inverse Quadratic
/// RBF.
///
/// @param x Numeric vector. The distances you wish to apply the Gaussian kernel
/// onto.
/// @param epsilon Float. Epsilon parameter for the RBF.
/// @param rbf_type String. Needs to be from `c("gaussian", "bump", "inverse_quadratic")`.
///
/// @return The affinities after the Kernel was applied.
///
/// @export
#[extendr]
fn rs_rbf_function(x: &[f64], epsilon: f64, rbf_type: &str) -> extendr_api::Result<Vec<f64>> {
    let rbf_fun = parse_rbf_types(rbf_type)
        .ok_or_else(|| extendr_api::Error::Other(format!("Invalid RBF function: {}", rbf_type)))?;

    let res: Vec<f64> = match rbf_fun {
        RbfType::Gaussian => rbf_gaussian(x, &epsilon),
        RbfType::Bump => rbf_bump(x, &epsilon),
        RbfType::InverseQuadratic => rbf_inverse_quadratic(x, &epsilon),
    };

    Ok(res)
}

/// Apply a Radial Basis Function (to a matrix)
///
/// @description Applies a radial basis function (RBF) to a given distance
/// matrix. Has at the option to apply a Gaussian, Bump or Inverse Quadratic
/// RBF.
///
/// @param x Numeric Matrix. The distances you wish to apply the Gaussian kernel
/// onto.
/// @param epsilon Float. Epsilon parameter for the RBF.
/// @param rbf_type String. Needs to be from `c("gaussian", "bump", "inverse_quadratic")`.
///
/// @return The affinities after the Kernel was applied.
///
/// @export
#[extendr]
fn rs_rbf_function_mat(
    x: RMatrix<f64>,
    epsilon: f64,
    rbf_type: &str,
) -> extendr_api::Result<extendr_api::RArray<f64, [usize; 2]>> {
    let x = r_matrix_to_faer(&x);

    let rbf_fun = parse_rbf_types(rbf_type)
        .ok_or_else(|| extendr_api::Error::Other(format!("Invalid RBF function: {}", rbf_type)))?;

    let res: Mat<f64> = match rbf_fun {
        RbfType::Gaussian => rbf_gaussian_mat(x, &epsilon),
        RbfType::Bump => rbf_bump_mat(x, &epsilon),
        RbfType::InverseQuadratic => rbf_inverse_quadratic_mat(x, &epsilon),
    };

    let res = faer_to_r_matrix(res.as_ref());

    Ok(res)
}

/// Calculates the TOM over an affinity matrix
///
/// @description Calculates the topological overlap measure for a given affinity matrix
/// x. Has the option to calculate the signed and unsigned version.
///
/// @param x Numerical matrix. Affinity matrix.
/// @param tom_type String. One of `c("v1", "v2")` - pending on choice, a different
/// normalisation method will be used.
/// @param signed Boolean. Shall the signed TOM be calculated. If set to `FALSE`, values
/// should be ≥ 0.
///
/// @return Returns the TOM matrix.
///
/// @export
#[extendr]
fn rs_tom(
    x: RMatrix<f64>,
    tom_type: &str,
    signed: bool,
) -> extendr_api::Result<extendr_api::RArray<f64, [usize; 2]>> {
    let x = r_matrix_to_faer(&x);
    let tom_version = parse_tom_types(tom_type).ok_or_else(|| {
        extendr_api::Error::Other(format!("Invalid TOM version type: {}", tom_type))
    })?;

    let tom_mat: Mat<f64> = match tom_version {
        TomType::Version1 => calc_tom(x, signed),
        TomType::Version2 => calc_tom_v2(x, signed),
    };

    let res = faer_to_r_matrix(tom_mat.as_ref());

    Ok(res)
}

/// Apply a range normalisation on a vector.
///
/// @description Applies a range normalisation on an R vector.
///
/// @param x Numerical vector. The data to normalise.
/// @param max_val Numeric. The upper bound value to normalise into. If set to 1,
/// the function will be equal to a min-max normalisation.
/// @param min_val Numeric. The lower bound value to normalise into. If set to 0,
/// the function will equal a min-max normalisation.
///
/// @return Normalised values
///
/// @export
#[extendr]
fn rs_range_norm(x: &[f64], max_val: f64, min_val: f64) -> Vec<f64> {
    let (x_min, x_max) = array_f64_max_min(x);
    let denom = x_max - x_min;
    let scale = (max_val - min_val) / denom;
    x.par_iter()
        .map(|x| (x - x_min) * scale + min_val)
        .collect()
}

extendr_module! {
    mod fun_helpers;
    fn rs_upper_triangle_to_dense;
    fn rs_coremo_quality;
    fn rs_ot_harmonic_sum;
    fn rs_rbf_function;
    fn rs_rbf_function_mat;
    fn rs_tom;
    fn rs_range_norm;
}
