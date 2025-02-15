use extendr_api::prelude::*;
use crate::utils_r_rust::r_matrix_to_faer;
use crate::helpers_linalg::column_covariance;


/// @export
#[extendr]
fn rs_covariance(
  x: RMatrix<f64>
) -> extendr_api::RArray<f64, [usize; 2]> {
  let mat = r_matrix_to_faer(x);

  let covar = column_covariance(&mat);

  let ncol = covar.ncols();

  RArray::new_matrix(ncol, ncol, |row, column| covar[(row, column)])
}

extendr_module! {
    mod fun_linalg;
    fn rs_covariance;
}