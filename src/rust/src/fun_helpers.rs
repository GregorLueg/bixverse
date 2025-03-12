use extendr_api::prelude::*;
use faer::Mat;
use crate::utils_r_rust::faer_to_r_matrix;


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
/// @return 
/// 
/// @export
#[extendr]
fn rs_upper_triangle_to_dense(
  cor_vector: Vec<f64>,
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

  faer_to_r_matrix(mat)
}

extendr_module! {
    mod fun_helpers;
    fn rs_upper_triangle_to_dense;
}