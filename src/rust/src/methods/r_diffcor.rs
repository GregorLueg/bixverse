use bixverse_rs::core::base::cors_similarity::column_pairwise_cor;
use bixverse_rs::methods::diffcor::*;
use bixverse_rs::prelude::*;
use extendr_api::prelude::*;

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_diffcor;
    fn rs_differential_cor;
}

///////////////
// Functions //
///////////////

/// Calculate the column wise differential correlation between two sets of data.
///
/// @description
/// `r lifecycle::badge("experimental")`
/// This function calculates the differential correlation based on the Fisher
/// method. For speed purposes, the function will only calculate the
/// differential correlation on the upper triangle of the two correlation
/// matrices.
///
/// @param x_a Numeric matrix a, samples x features.
/// @param x_b Numeric matrix b, samples x features. Needs the same number of
/// columns as `x_a`.
/// @param spearman Boolean. Shall the Spearman correlation be calculated
/// instead of Pearson.
///
/// @returns A list containing, one entry per upper-triangle feature pair:
///  \itemize{
///   \item r_a - The correlation coefficients in the upper triangle of
///   matrix a.
///   \item r_b - The correlation coefficients in the upper triangle of
///   matrix b.
///   \item z_score - The z-scores of the difference in correlation
///   coefficients.
///   \item p_val - The z-scores transformed to two-sided p-values.
/// }
///
/// @export
#[extendr]
fn rs_differential_cor(
    x_a: RMatrix<f64>,
    x_b: RMatrix<f64>,
    spearman: bool,
) -> extendr_api::Result<List> {
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

    let cor_a = column_pairwise_cor(&mat_a, spearman);
    let cor_b = column_pairwise_cor(&mat_b, spearman);

    let diff_cor =
        calculate_diff_correlation::<f64>(&cor_a, &cor_b, n_sample_a, n_sample_b, spearman);

    Ok(list!(
        r_a = diff_cor.r_a,
        r_b = diff_cor.r_b,
        z_score = diff_cor.z_score,
        p_val = diff_cor.p_vals
    ))
}
