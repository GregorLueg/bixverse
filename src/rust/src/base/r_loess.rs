use extendr_api::*;

use bixverse_rs::core::base::loess::*;

////////////////////
// extendr Module //
////////////////////

extendr_module! {
  mod r_loess;
  fn rs_2d_loess;
}

///////////////
// Functions //
///////////////

/// Rust implementation of a Loess function
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Fits a Loess regression of `y` on `x`. Points where either value is
/// non-finite are dropped from the fit.
///
/// @param x Numeric. The x values to fit.
/// @param y Numeric. The y values to fit.
/// @param span Numeric. The span parameter. Needs to be in `(0, 1]`.
/// @param degree Integer. Either 1 (linear) or 2 (quadratic). Other values
/// will cause an error.
///
/// @returns A list with the following items
/// \itemize{
///   \item predicted - The predicted values, `0` for dropped points.
///   \item residuals - The residuals for every data point, `0` for dropped
///   points.
///   \item valid_idx - 1-based indices of the points included in the fit.
/// }
///
/// @export
#[extendr]
fn rs_2d_loess(x: &[f64], y: &[f64], span: f64, degree: usize) -> List {
    let loess = LoessRegression::new(span, degree);
    let loess_res: LoessRes<f64> = loess.fit(x, y);

    list!(
        predicted = loess_res.fitted_vals,
        residuals = loess_res.residuals,
        valid_idx = loess_res
            .valid_indices
            .iter()
            .map(|x| (*x + 1) as i32)
            .collect::<Vec<i32>>(),
    )
}
