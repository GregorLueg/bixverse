use crate::assert_same_len;

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

        let mut sorted_points = valid.clone();
        sorted_points.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        let x_sorted: Vec<f64> = sorted_points.iter().map(|(_, x, _)| *x).collect();
        let y_sorted: Vec<f64> = sorted_points.iter().map(|(_, _, y)| *y).collect();

        let n_valid = valid.len();
        let no_neighbours = ((n_valid as f64) * self.span).max(1.0) as usize;

        let mut fitted_values = vec![0.0; n];
        let mut residuals = vec![0.0; n];

        for (orig_idx, x_val, y_val) in &valid {
            let fitted_val = self.fit_local_regression(x_val, &x_sorted, &y_sorted, no_neighbours);

            fitted_values[*orig_idx] = fitted_val;
            residuals[*orig_idx] = y_val - fitted_val;
        }

        LoessRes {
            fitted_vals: fitted_values,
            residuals,
            valid_indices: valid.iter().map(|(idx, _, _)| *idx).collect(),
        }
    }

    /// Helper to fit the local regression
    ///
    /// ### Params
    ///
    /// * `target_x` - The target variable at this point
    /// * `x` - Slice of the predictor variable
    /// * `y` - Slice of the response variable
    /// * `k` - Number of neighbours to use
    ///
    /// ### Returns
    ///
    /// The value for the target x
    fn fit_local_regression(&self, target_x: &f64, x: &[f64], y: &[f64], k: usize) -> f64 {
        let neighbours = self.find_neighbors(target_x, x, k);

        // early return if there are no neighbours
        if neighbours.is_empty() {
            return 0.0;
        }

        let max_dist = neighbours
            .iter()
            .map(|&i| (x[i] - target_x).abs())
            .fold(0.0, f64::max);

        // all neighbors are at the same x value
        if max_dist == 0.0 {
            let weights_sum: f64 = neighbours.len() as f64;
            return neighbours.iter().map(|&i| y[i]).sum::<f64>() / weights_sum;
        }

        let weights: Vec<f64> = neighbours
            .iter()
            .map(|&i| self.tricube_weight((x[i] - target_x).abs() / max_dist))
            .collect();

        self.weighted_polynomial_fit(
            target_x,
            &neighbours.iter().map(|&i| x[i]).collect::<Vec<_>>(),
            &neighbours.iter().map(|&i| y[i]).collect::<Vec<_>>(),
            &weights,
        )
    }

    /// Helper function get k nearest neighbours
    ///
    /// ### Params
    ///
    /// * `target_x` - Target value.
    /// * `x` - Remaining other values.
    /// * `k` - Number of neighbours.
    ///
    /// ### Returns
    ///
    /// The neighbour positions
    fn find_neighbors(&self, target_x: &f64, x: &[f64], k: usize) -> Vec<usize> {
        let mut distances: Vec<(usize, f64)> = x
            .iter()
            .enumerate()
            .map(|(i, x)| (i, (x - target_x).abs()))
            .collect();

        distances.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        distances.truncate(k);

        distances.into_iter().map(|(i, _)| i).collect()
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
    fn tricube_weight(&self, u: f64) -> f64 {
        if u.abs() >= 1.0 {
            0.0
        } else {
            let temp = 1.0 - u.abs().powi(3);
            temp.powi(3)
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
    fn weighted_polynomial_fit(
        &self,
        target_x: &f64,
        x_vals: &[f64],
        y_vals: &[f64],
        weights: &[f64],
    ) -> f64 {
        let n = x_vals.len();
        if n == 0 {
            return 0.0;
        }

        match self.loess_type {
            LoessFunc::Linear => self.weighted_linear_fit(target_x, x_vals, y_vals, weights),
            LoessFunc::Quadratic => self.weighted_quadratic_fit(target_x, x_vals, y_vals, weights),
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
        let w_sum = w.iter().sum::<f64>();
        if w_sum == 0_f64 {
            let res = y.iter().sum::<f64>() / y.len() as f64;
            return res;
        }

        let wx_sum = w.iter().zip(x).map(|(w, x)| w * x).sum::<f64>();
        let wy_sum = w.iter().zip(y).map(|(w, y)| w * y).sum::<f64>();
        let wxx_sum = w.iter().zip(x).map(|(w, x)| w * x * x).sum::<f64>();
        let wxy_sum = w
            .iter()
            .zip(x.iter().zip(y.iter()))
            .map(|(w, (x, y))| w * x * y)
            .sum::<f64>();

        let x_mean = wx_sum / w_sum;
        let y_mean = wy_sum / w_sum;

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
        let n = x.len();
        // Fall back to linear for insufficient points
        if n < 3 {
            return self.weighted_linear_fit(target_x, x, y, w);
        }

        let mut a = [[0.0; 3]; 3]; // 3x3 matrix
        let mut b = [0.0; 3];

        for i in 0..n {
            let x_i = x[i];
            let y_i = y[i];
            let w_i = w[i];

            let basis = [1.0, x_i, x_i * x_i];

            for j in 0..3 {
                for k in 0..3 {
                    a[j][k] += w_i * basis[j] * basis[k];
                }
                b[j] += w_i * basis[j] * y_i;
            }
        }

        match solve_3x3_system(&a, &b) {
            Some(coeffs) => {
                // Evaluate polynomial at target_x
                coeffs[0] + coeffs[1] * target_x + coeffs[2] * target_x * target_x
            }
            None => {
                // Fall back to linear regression
                self.weighted_linear_fit(target_x, x, y, w)
            }
        }
    }
}
