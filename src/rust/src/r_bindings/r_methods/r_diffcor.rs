use extendr_api::prelude::*;

use crate::core::base::cors_similarity::column_cor;
use crate::core::methods::diffcor::*;
use crate::utils::r_rust_interface::*;

/// Calculate the column wise differential correlation between two sets of data.
///
/// @description This function calculates the differential correlation based on
/// the Fisher method. For speed purposes, the function will only calculate the
/// differential correlation on the upper triangle of the two correlation
/// matrices.
/// WARNING! Incorrect use can cause kernel crashes. Wrapper around the Rust
/// functions with type checks are provided in the package.
///
/// @param x_a R matrix a to be used for the differential correlation analysis.
/// @param x_b R matrix a to be used for the differential correlation analysis.
/// @param spearman Shall the Spearman correlation be calculated instead of
/// Pearson.
///
/// @return A list containing:
///  \itemize{
///   \item r_a - The correlation coefficients in the upper triangle of
///   matrix a.
///   \item r_b - The correlation coefficients in the upper triangle of
///   matrix b.
///   \item z_score - The z-scores of the difference in correlation
///   coefficients.
///   \item p_val - The z-scores transformed to p-values.
/// }
///
/// @export
#[extendr]
fn rs_differential_cor(x_a: RMatrix<f64>, x_b: RMatrix<f64>, spearman: bool) -> List {
    assert!(
        x_a.ncols() == x_b.ncols(),
        "Input matrices must have the same number of columns. Found {} columns in first matrix and {} in second.",
        x_a.ncols(),
        x_b.ncols(),
      );
    let n_sample_a = x_a.nrows();
    let n_sample_b = x_b.nrows();
    let mat_a = r_matrix_to_faer(&x_a);
    let mat_b = r_matrix_to_faer(&x_b);

    let cor_a = column_cor(&mat_a, spearman);
    let cor_b = column_cor(&mat_b, spearman);

    let diff_cor = calculate_diff_correlation(&cor_a, &cor_b, n_sample_a, n_sample_b, spearman);

    list!(
        r_a = diff_cor.r_a,
        r_b = diff_cor.r_b,
        z_score = diff_cor.z_score,
        p_val = diff_cor.p_vals
    )
}

extendr_module! {
    mod r_diffcor;
    fn rs_differential_cor;
}
