use extendr_api::*;

use crate::core::base::loess::*;

/// Rust implementation of a Loess function
///
/// @param x Numeric. The x values to fit.
/// @param y Numeric. The y values to fit.
/// @param span Numeric. The span parameter. Needs to be between 0.1 and 1.
/// @param degree Integer. Either 1 (linear) or 2 (quadratic). Other values
/// will cause an error.
///
/// @return A list with the following items
/// \itemize{
///   \item predicted - The predicted values.
///   \item residuals - The residuals for every data point.
///   \item valid_idx - Which data indices were included.
/// }
///
/// @export
#[extendr]
fn rs_2d_loess(x: &[f64], y: &[f64], span: f64, degree: usize) -> List {
    let loess = LoessRegression::new(span, degree);
    let loess_res: LoessRes = loess.fit(x, y);

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

extendr_module! {
  mod r_loess;
  fn rs_2d_loess;
}
