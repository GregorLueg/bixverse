use extendr_api::prelude::*;
use faer::Mat;
use rayon::prelude::*;

use crate::core::base::stats::calculate_critval;
use crate::utils::general::array_max_min;
use crate::utils::r_rust_interface::*;

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

/// Generate a vector-based representation of the upper triangle of a matrix
///
/// @description This function generates a vector from the upper triangle of
/// a given symmetric matrix. You have the option to remove the diagonal with
/// setting shift to 1.
///
/// @param x Numeric vector. The vector of correlation coefficients that you
/// want to use to go back to a dense matrix.
/// @param shift Integer. If you want to apply a shift, i.e. included the diagonal
/// values = 0; or excluded the diagonal values = 1.
///
/// @return The dense R matrix.
///
/// @export
#[extendr]
fn rs_dense_to_upper_triangle(x: RMatrix<f64>, shift: usize) -> Vec<f64> {
    let n = x.ncols();

    let total_elements = if shift == 0 {
        n * (n + 1) / 2
    } else {
        n * (n - 1) / 2
    };

    let mut vals = Vec::with_capacity(total_elements);

    for i in 0..n {
        let start_j = i + shift;
        for j in start_j..n {
            vals.push(x[[i, j]]);
        }
    }

    vals
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
    let (x_min, x_max) = array_max_min(x);
    let denom = x_max - x_min;
    let scale = (max_val - min_val) / denom;
    x.par_iter()
        .map(|x| (x - x_min) * scale + min_val)
        .collect()
}

/// Calculate the critical value
///
/// This function calculates the critical value for a given set based on random
/// permutations and a given alpha value.
///
/// @param values Numeric vector. The full data set for which to calculate the
/// critical value.
/// @param iters Integer. Number of random permutations to use.
/// @param alpha Float. The alpha value. For example, 0.001 would mean that the
/// critical value is smaller than 0.1 percentile of the random permutations.
/// @param seed Integer. For reproducibility purposes
///
/// @return The critical value for the given parameters.
///
/// @export
#[extendr]
fn rs_critval(values: &[f64], iters: usize, alpha: f64, seed: usize) -> f64 {
    calculate_critval(values, iters, &alpha, seed)
}

/// Calculate the critical value
///
/// This function calculates the critical value for a given set based on random
/// permutations and a given alpha value.
///
/// @param mat Numeric matrix. The (symmetric matrix with all of the values).
/// @param iters Integer. Number of random permutations to use.
/// @param alpha Float. The alpha value. For example, 0.001 would mean that the
/// critical value is smaller than 0.1 percentile of the random permutations.
/// @param seed Integer. For reproducibility purposes
///
/// @return The critical value for the given parameters.
///
/// @export
#[extendr]
fn rs_critval_mat(mat: RMatrix<f64>, iters: usize, alpha: f64, seed: usize) -> f64 {
    let mut values: Vec<f64> = Vec::new();
    for r in 0..mat.nrows() {
        for c in (r + 1)..mat.ncols() {
            values.push(mat[[r, c]]);
        }
    }

    calculate_critval(&values, iters, &alpha, seed)
}

extendr_module! {
    mod r_helpers;
    fn rs_upper_triangle_to_dense;
    fn rs_dense_to_upper_triangle;
    fn rs_ot_harmonic_sum;
    fn rs_range_norm;
    fn rs_critval;
    fn rs_critval_mat;
}
