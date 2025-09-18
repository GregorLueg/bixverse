use crate::assert_same_len;
use rayon::prelude::*;

/// Structure to store the Loess results
///
/// ### Params
///
/// * `fitted_vals` - The values fitted by the function.
/// * `residuals` - The residuals.
/// * `valid_indices` - Which index positions were valid.
#[derive(Debug, Clone)]
pub struct LoessRes {
    pub fitted_vals: Vec<f64>,
    pub residuals: Vec<f64>,
    pub valid_indices: Vec<usize>,
}

#[derive(Debug, Clone)]
pub enum LoessFunc {
    /// Linear version of the Loess function
    Linear,
    /// Quadratic version of the Loess function
    Quadratic,
}

/// Parse the type of Loess function
///
/// ### Params
///
/// * `option` - Usize defining the degrees of freedom
///
/// ### Return
///
/// The option of the `LoessFunc`
pub fn parse_loess_fun(option: &usize) -> Option<LoessFunc> {
    match option {
        1 => Some(LoessFunc::Linear),
        2 => Some(LoessFunc::Quadratic),
        _ => None,
    }
}

/////////////
// Helpers //
/////////////

fn solve_3x3_system(a: &[[f64; 3]; 3], b: &[f64; 3]) -> Option<[f64; 3]> {
    let mut matrix = *a;
    let mut rhs = *b;

    for i in 0..3 {
        // Find pivot
        let mut pivot_row = i;
        for j in (i + 1)..3 {
            if matrix[j][i].abs() > matrix[pivot_row][i].abs() {
                pivot_row = j;
            }
        }

        // Swap rows if needed
        if pivot_row != i {
            matrix.swap(i, pivot_row);
            rhs.swap(i, pivot_row);
        }

        // Check for singular matrix
        if matrix[i][i].abs() < 1e-12 {
            return None;
        }

        // Eliminate
        for j in (i + 1)..3 {
            let factor = matrix[j][i] / matrix[i][i];
            for k in i..3 {
                matrix[j][k] -= factor * matrix[i][k];
            }
            rhs[j] -= factor * rhs[i];
        }
    }

    let mut solution = [0.0; 3];
    for i in (0..3).rev() {
        solution[i] = rhs[i];
        for j in (i + 1)..3 {
            solution[i] -= matrix[i][j] * solution[j];
        }
        solution[i] /= matrix[i][i];
    }

    Some(solution)
}

///////////
// Loess //
///////////

pub struct LoessRegression {
    span: f64,
    loess_type: LoessFunc,
}

impl LoessRegression {
    /// Generate a new instance of the Loess regression
    ///
    /// ### Params
    ///
    /// * `span` -
    /// * `degree` -
    ///
    /// ### Return
    ///
    /// Initialised class
    pub fn new(span: f64, degree: usize) -> Self {
        assert!(span > 0.0 && span <= 1.0, "Span must be between 0 and 1");
        assert!(
            degree == 1 || degree == 2,
            "Only linear (1) and quadratic (2) supported"
        );

        let loess_type: LoessFunc = parse_loess_fun(&degree).unwrap();

        Self { span, loess_type }
    }

    /// Fit the loess function (for a two variable system)
    ///
    /// ### Params
    ///
    /// * `x` - The response variable
    /// * `y` - The predictor variable
    ///
    /// ### Returns
    ///
    /// The fit results in form of a `LoessRes`
    pub fn fit<T>(&self, x: &[T], y: &[T]) -> LoessRes
    where
        T: Copy + Into<f64> + PartialOrd,
        f64: From<T>,
    {
        assert_same_len!(x, y);

        let n = x.len();
        let valid: Vec<(usize, f64, f64)> = x
            .iter()
            .zip(y.iter())
            .enumerate()
            .filter_map(|(i, (&x, &y))| {
                let x_f64: f64 = x.into();
                let y_f64: f64 = y.into();
                if x_f64.is_finite() && y_f64.is_finite() {
                    Some((i, x_f64, y_f64))
                } else {
                    None
                }
            })
            .collect();

        if valid.is_empty() {
            return LoessRes {
                fitted_vals: vec![0.0; n],
                residuals: vec![0.0; n],
                valid_indices: Vec::new(),
            };
        }

        let mut sorted_points = valid;
        sorted_points.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        let n_valid = sorted_points.len();
        let no_neighbours = ((n_valid as f64) * self.span).max(1.0) as usize;

        let mut fitted_values = vec![0.0; n];
        let mut residuals = vec![0.0; n];

        // parallelise for speed up
        let results: Vec<_> = sorted_points
            .par_iter()
            .map(|(orig_idx, x_val, y_val)| {
                let fitted_val = self.fit_point(&sorted_points, *x_val, no_neighbours);
                (*orig_idx, fitted_val, y_val - fitted_val)
            })
            .collect();

        for (orig_idx, fitted_val, residual) in results {
            fitted_values[orig_idx] = fitted_val;
            residuals[orig_idx] = residual;
        }

        LoessRes {
            fitted_vals: fitted_values,
            residuals,
            valid_indices: sorted_points.iter().map(|(idx, _, _)| *idx).collect(),
        }
    }

    /// Fits a given point
    ///
    /// ### Params
    ///
    /// * `sorted_points` - A slice of tuples of the position, x and y value
    /// * `target_x` - The target value
    /// * `k` - Number of neighbours
    ///
    /// ### Returns
    ///
    /// The fitted value
    fn fit_point(&self, sorted_points: &[(usize, f64, f64)], target_x: f64, k: usize) -> f64 {
        let neighbors = self.find_neighbors_binary(sorted_points, target_x, k);

        if neighbors.is_empty() {
            return 0.0;
        }

        let max_dist = neighbors
            .iter()
            .map(|&i| (sorted_points[i].1 - target_x).abs())
            .fold(0.0, f64::max);

        if max_dist == 0.0 {
            return neighbors.iter().map(|&i| sorted_points[i].2).sum::<f64>()
                / neighbors.len() as f64;
        }

        let mut x_vals = [0.0; 64];
        let mut y_vals = [0.0; 64];
        let mut weights = [0.0; 64];

        let n_neighbors = neighbors.len().min(64);
        let inv_max_dist = 1.0 / max_dist;

        for (i, &idx) in neighbors.iter().take(n_neighbors).enumerate() {
            let (_, nx, ny) = sorted_points[idx];
            x_vals[i] = nx;
            y_vals[i] = ny;
            weights[i] = self.tricube_weight((nx - target_x).abs() * inv_max_dist);
        }

        self.weighted_polynomial_fit(
            &target_x,
            &x_vals[..n_neighbors],
            &y_vals[..n_neighbors],
            &weights[..n_neighbors],
        )
    }

    /// Find neighhbours via binary search
    ///
    /// Uses binary search under the hood for speed.
    ///
    /// ### Params
    ///
    /// `sorted_points` - A slice of tuples of the position, x and y value
    /// `target_x` - The target value
    /// `k` - Number of neighbour
    fn find_neighbors_binary(
        &self,
        sorted_points: &[(usize, f64, f64)],
        target_x: f64,
        k: usize,
    ) -> Vec<usize> {
        let n = sorted_points.len();
        if k >= n {
            return (0..n).collect();
        }

        // Binary search for insertion point
        let insert_pos = sorted_points
            .binary_search_by(|probe| probe.1.partial_cmp(&target_x).unwrap())
            .unwrap_or_else(|pos| pos);

        // Expand around insertion point
        let mut l = insert_pos;
        let mut r = insert_pos;
        let mut neighbors = Vec::with_capacity(k);

        for _ in 0..k {
            let left_dist = if l > 0 {
                (sorted_points[l - 1].1 - target_x).abs()
            } else {
                f64::INFINITY
            };
            let right_dist = if r < n {
                (sorted_points[r].1 - target_x).abs()
            } else {
                f64::INFINITY
            };

            if left_dist <= right_dist && l > 0 {
                l -= 1;
                neighbors.push(l);
            } else if r < n {
                neighbors.push(r);
                r += 1;
            } else {
                break;
            }
        }

        neighbors
    }

    /// Tricube weight function: (1 - |u|³)³ for |u| < 1, 0 otherwise
    ///
    /// ### Params
    ///
    /// * `u` - The value
    ///
    /// ### Returns
    ///
    /// The tricube weight
    #[inline]
    fn tricube_weight(&self, u: f64) -> f64 {
        if u >= 1.0 {
            0.0
        } else {
            let temp = 1.0 - u * u * u;
            temp * temp * temp
        }
    }

    /// Weighted fit helper
    ///
    /// ### Params
    ///
    /// * `target_x` - The target variable at this point
    /// * `x` - Slice of the predictor variable
    /// * `y` - Slice of the response variable
    /// * `w` - Slice of the weights
    ///
    /// ### Returns
    ///
    /// Value at this position
    fn weighted_polynomial_fit(&self, target_x: &f64, x: &[f64], y: &[f64], w: &[f64]) -> f64 {
        match self.loess_type {
            LoessFunc::Linear => self.weighted_linear_fit(target_x, x, y, w),
            LoessFunc::Quadratic => self.weighted_quadratic_fit(target_x, x, y, w),
        }
    }

    /// Helper function for linear fits
    ///
    /// ### Params
    ///
    /// * `target_x` - The target variable at this point
    /// * `x` - Slice of the predictor variable
    /// * `y` - Slice of the response variable
    /// * `w` - Slice of the weights
    ///
    /// ### Return
    ///
    /// Value at this position with a linear regression
    fn weighted_linear_fit(&self, target_x: &f64, x: &[f64], y: &[f64], w: &[f64]) -> f64 {
        let mut w_sum = 0.0;
        let mut wx_sum = 0.0;
        let mut wy_sum = 0.0;
        let mut wxx_sum = 0.0;
        let mut wxy_sum = 0.0;

        for i in 0..x.len() {
            let wi = w[i];
            let xi = x[i];
            let yi = y[i];

            w_sum += wi;
            wx_sum += wi * xi;
            wy_sum += wi * yi;
            wxx_sum += wi * xi * xi;
            wxy_sum += wi * xi * yi;
        }

        if w_sum == 0.0 {
            return y.iter().sum::<f64>() / y.len() as f64;
        }

        let inv_w_sum = 1.0 / w_sum;
        let x_mean = wx_sum * inv_w_sum;
        let y_mean = wy_sum * inv_w_sum;

        let numerator = wxy_sum - w_sum * x_mean * y_mean;
        let denominator = wxx_sum - w_sum * x_mean * x_mean;

        if denominator.abs() < 1e-12 {
            return y_mean;
        }

        let slope = numerator / denominator;
        let intercept = y_mean - slope * x_mean;

        intercept + slope * target_x
    }

    /// Helper function for quadratic fits
    ///
    /// ### Params
    ///
    /// * `target_x` - The target variable at this point
    /// * `x` - Slice of the predictor variable
    /// * `y` - Slice of the response variable
    /// * `w` - Slice of the weights
    ///
    /// ### Return
    ///
    /// Value at this position with a quadratic regression
    fn weighted_quadratic_fit(&self, target_x: &f64, x: &[f64], y: &[f64], w: &[f64]) -> f64 {
        if x.len() < 3 {
            return self.weighted_linear_fit(target_x, x, y, w);
        }

        let mut a = [[0.0; 3]; 3];
        let mut b = [0.0; 3];

        for i in 0..x.len() {
            let xi = x[i];
            let yi = y[i];
            let wi = w[i];
            let xi2 = xi * xi;

            // Manual matrix construction (faster than loops)
            a[0][0] += wi;
            a[0][1] += wi * xi;
            a[0][2] += wi * xi2;
            a[1][1] += wi * xi2;
            a[1][2] += wi * xi * xi2;
            a[2][2] += wi * xi2 * xi2;

            b[0] += wi * yi;
            b[1] += wi * xi * yi;
            b[2] += wi * xi2 * yi;
        }

        // Symmetric matrix
        a[1][0] = a[0][1];
        a[2][0] = a[0][2];
        a[2][1] = a[1][2];

        match solve_3x3_system(&a, &b) {
            Some(coeffs) => coeffs[0] + coeffs[1] * target_x + coeffs[2] * target_x * target_x,
            None => self.weighted_linear_fit(target_x, x, y, w),
        }
    }
}
