use num_traits::Float;
use std::iter::Sum;

/// Simple linear regression
///
/// Fits y = b0 + b1 * x using ordinary least squares.
///
/// ### Params
///
/// * `x` - Independent variable
/// * `y` - Dependent variable
///
/// ### Returns
///
/// Tuple of (intercept, slope)
pub fn linear_regression<F>(x: &[F], y: &[F]) -> (F, F)
where
    F: Float + Sum,
{
    let n = F::from(x.len()).unwrap();
    let sum_x: F = x.iter().cloned().sum();
    let sum_y: F = y.iter().cloned().sum();
    let sum_xy: F = x.iter().zip(y).map(|(&xi, &yi)| xi * yi).sum();
    let sum_xx: F = x.iter().map(|&xi| xi * xi).sum();

    let slope = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x * sum_x);
    let intercept = (sum_y - slope * sum_x) / n;

    (intercept, slope)
}
