use faer::{
    linalg::solvers::{PartialPivLu, Solve},
    traits::AddByRef,
    ColRef, Mat, MatRef,
};
use rayon::iter::*;

///////////
// Enums //
///////////

/// Binning strategy enum
#[derive(Debug, Clone)]
pub enum BinningStrategy {
    /// Equal width bins (equal distance between bin edges)
    EqualWidth,
    /// Equal frequency bins (approximately equal number of values per bin)
    EqualFrequency,
}

////////////
// Params //
////////////

/// Parsing the binning strategy
///
/// ### Params
///
/// * `s` - string defining the binning strategy
///
/// ### Returns
///
/// The `BinningStrategy`.
pub fn parse_bin_strategy_type(s: &str) -> Option<BinningStrategy> {
    match s.to_lowercase().as_str() {
        "equal_width" => Some(BinningStrategy::EqualWidth),
        "equal_freq" => Some(BinningStrategy::EqualFrequency),
        _ => None,
    }
}

////////////////////
// Util functions //
////////////////////

/// Generate the rank of a vector with tie correction.
///
/// ### Params
///
/// * `vec` - The slice of numericals to rank.
///
/// ### Returns
///
/// The ranked vector (also f64)
pub fn rank_vector<T>(vec: &[T]) -> Vec<f64>
where
    T: Copy + PartialOrd + PartialEq,
{
    let n = vec.len();
    if n == 0 {
        return Vec::new();
    }

    let mut indexed_values: Vec<(T, usize)> = vec
        .iter()
        .copied()
        .enumerate()
        .map(|(i, v)| (v, i))
        .collect();

    indexed_values
        .sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let mut ranks = vec![0.0; n];
    let mut i = 0;
    while i < n {
        let current_value = indexed_values[i].0;
        let start = i;
        while i < n && indexed_values[i].0 == current_value {
            i += 1;
        }
        let avg_rank = (start + i + 1) as f64 / 2.0;
        for j in start..i {
            ranks[indexed_values[j].1] = avg_rank;
        }
    }
    ranks
}

/// Calculates the columns means of a matrix
///
/// ### Params
///
/// * `mat` - The matrix for which to calculate the column-wise means
///
/// ### Returns
///
/// Vector of the column means.
pub fn col_means(mat: MatRef<f64>) -> Vec<f64> {
    let n_rows = mat.nrows();
    let ones = Mat::from_fn(n_rows, 1, |_, _| 1.0);
    let means = (ones.transpose() * mat) / n_rows as f64;

    means.row(0).iter().cloned().collect()
}

/// Calculates the column sums of a matrix
///
/// ### Params
///
/// * `mat` - The matrix for which to calculate the column-wise sums
///
/// ### Returns
///
/// Vector of the column sums.
pub fn col_sums(mat: MatRef<f64>) -> Vec<f64> {
    let n_rows = mat.nrows();
    let ones = Mat::from_fn(n_rows, 1, |_, _| 1.0);
    let col_sums = ones.transpose() * mat;

    col_sums.row(0).iter().cloned().collect()
}

/// Calculate the column standard deviations
///
/// ### Params
///
/// * `mat` - The matrix for which to calculate the column-wise standard deviations
///
/// ### Returns
///
/// Vector of the column standard deviations.
pub fn col_sds(mat: MatRef<f64>) -> Vec<f64> {
    let n = mat.nrows() as f64;
    let n_cols = mat.ncols();

    // Calculate means and SDs in one pass
    let (_, m2): (Vec<f64>, Vec<f64>) = (0..n_cols)
        .map(|j| {
            let mut mean = 0.0;
            let mut m2 = 0.0;
            let mut count = 0.0;

            for i in 0..mat.nrows() {
                count += 1.0;
                let delta = mat[(i, j)] - mean;
                mean += delta / count;
                let delta2 = mat[(i, j)] - mean;
                m2 += delta * delta2;
            }
            (mean, (m2 / (n - 1.0)).sqrt())
        })
        .unzip();

    m2
}

/// Scale a matrix
///
/// ### Params
///
/// * `mat` - The matrix on which to apply column-wise scaling
/// * `scale_sd` - Shall the standard deviation be equalised across columns
///
/// ### Returns
///
/// The scaled matrix.
pub fn scale_matrix_col(mat: &MatRef<f64>, scale_sd: bool) -> Mat<f64> {
    let n_rows = mat.nrows();
    let n_cols = mat.ncols();

    let mut means = vec![0.0; n_cols];
    for j in 0..n_cols {
        for i in 0..n_rows {
            means[j] += mat[(i, j)];
        }
        means[j] /= n_rows as f64;
    }

    let mut result = mat.to_owned();
    for j in 0..n_cols {
        let mean = means[j];
        for i in 0..n_rows {
            result[(i, j)] -= mean;
        }
    }

    if !scale_sd {
        return result;
    }

    let mut std_devs = vec![0.0; n_cols];
    for j in 0..n_cols {
        for i in 0..n_rows {
            let val = result[(i, j)];
            std_devs[j] += val * val;
        }
        std_devs[j] = (std_devs[j] / (n_rows as f64 - 1.0)).sqrt();
        if std_devs[j] < 1e-10 {
            std_devs[j] = 1.0;
        }
    }

    for j in 0..n_cols {
        let std_dev = std_devs[j];
        for i in 0..n_rows {
            result[(i, j)] /= std_dev;
        }
    }

    result
}

/// Column wise L2 normalisation
///
/// ### Params
///
/// * `mat` - The matrix on which to apply column-wise L2 normalisation
///
/// ### Returns
///
/// The matrix with the columns being L2 normalised.
pub fn normalise_matrix_col_l2(mat: &MatRef<f64>) -> Mat<f64> {
    let mut normalized = mat.to_owned();

    for j in 0..mat.ncols() {
        let col = mat.col(j);
        let norm = col.norm_l2();

        if norm > 1e-10 {
            for i in 0..mat.nrows() {
                normalized[(i, j)] = mat[(i, j)] / norm;
            }
        }
    }

    normalized
}

/// Column wise rank normalisation
///
/// ### Params
///
/// * `mat` - The matrix on which to apply column-wise rank normalisation
///
/// ### Returns
///
/// The matrix with the columns being rank normalised.
pub fn rank_matrix_col(mat: &MatRef<f64>) -> Mat<f64> {
    let mut ranked_mat = Mat::zeros(mat.nrows(), mat.ncols());

    // Parallel ranking directly into the matrix
    ranked_mat
        .par_col_iter_mut()
        .enumerate()
        .for_each(|(col_idx, mut col)| {
            let original_col: Vec<f64> = mat.col(col_idx).iter().copied().collect();
            let ranks = rank_vector(&original_col);

            // Write ranks directly to the matrix column
            for (row_idx, &rank) in ranks.iter().enumerate() {
                col[row_idx] = rank;
            }
        });

    ranked_mat
}

/// Equal width binning for a single column
///
/// ### Params
///
/// * `col` - Column refence to bin with equal width strategy
/// * `n_bins` - Number of bins to use
///
/// ### Returns
///
/// Binned vector
fn bin_equal_width(col: &ColRef<f64>, n_bins: usize) -> Vec<usize> {
    let (min_val, max_val) = col
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(min, max), x| {
            (min.min(*x), max.max(*x))
        });

    let range = max_val - min_val;

    if range == 0.0 {
        return vec![0; col.nrows()];
    }

    let step = range / n_bins as f64;

    col.iter()
        .map(|x| {
            let bin = ((*x - min_val) / step).floor() as usize;
            bin.min(n_bins - 1)
        })
        .collect()
}

/// Equal frequency binning for a single column
///
/// ### Params
///
/// * `col` - Column refence to bin with equal frequency strategy
/// * `n_bins` - Number of bins to use
///
/// ### Returns
///
/// Binned vector
fn bin_equal_frequency(col: &faer::col::ColRef<f64>, n_bins: usize) -> Vec<usize> {
    let n_rows = col.nrows();

    // Create a copy of the column values for sorting
    let mut sorted_values: Vec<f64> = col.iter().copied().collect();
    sorted_values.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    // Calculate split points (quantiles) like infotheo does
    let mut split_points = Vec::with_capacity(n_bins);
    for k in 1..n_bins {
        let idx = (k * n_rows) / n_bins;
        // Use the value at the calculated index as split point
        split_points.push(sorted_values[idx.min(n_rows - 1)]);
    }

    // Assign bins based on split points
    let mut bin_assignments = vec![0; n_rows];
    for (i, &value) in col.iter().enumerate() {
        let mut bin = 0;
        for &split in &split_points {
            if value >= split {
                bin += 1;
            } else {
                break;
            }
        }
        bin_assignments[i] = bin.min(n_bins - 1);
    }

    bin_assignments
}

/// Column wise binning
///
/// ### Params
///
/// * `mat` - The matrix on which to apply column-wise binning
/// * `n_bins` - Optional number of bins. If not provided, will default to
///   `sqrt(nrow)`
///
/// ### Returns
///
/// The matrix with the columns being binned into equal distances.
pub fn bin_matrix_cols(
    mat: &MatRef<f64>,
    n_bins: Option<usize>,
    strategy: BinningStrategy,
) -> Mat<usize> {
    let (n_rows, n_cols) = mat.shape();
    let n_bins = n_bins.unwrap_or_else(|| (n_rows as f64).sqrt() as usize);

    let binned_vals: Vec<Vec<usize>> = mat
        .par_col_iter()
        .map(|col| match strategy {
            BinningStrategy::EqualWidth => bin_equal_width(&col, n_bins),
            BinningStrategy::EqualFrequency => bin_equal_frequency(&col, n_bins),
        })
        .collect();

    Mat::from_fn(n_rows, n_cols, |i, j| binned_vals[j][i])
}

////////////////////
// Matrix solvers //
////////////////////

/// Sylvester solver for three matrix systems
///
/// Solves a system of `AX + XB = C`. Pending on the size of the underlying
/// matrices, the algorithm will solve this directly or iteratively.
///
/// ### Params
///
/// * `mat_a` - Matrix A of the system
/// * `mat_b` - Matrix B of the system
/// * `mat_c` - Matrix C of the system
///
/// ### Returns
///
/// The matrix X
pub fn sylvester_solver(mat_a: &MatRef<f64>, mat_b: &MatRef<f64>, mat_c: &MatRef<f64>) -> Mat<f64> {
    let m = mat_a.nrows();
    let n = mat_b.ncols();

    if m * n < 1000 {
        // For small problems, use direct method
        sylvester_solver_direct(mat_a, mat_b, mat_c)
    } else {
        // For large problems use the iterative method
        sylvester_solver_iterative(mat_a, mat_b, mat_c, 50, 1e-6)
    }
}

/// Iterative Sylvester solver using fixed-point iteration
///
/// Solves a system of `AX + XB = C`. Uses an iterative approach more
/// appropriate for large matrix systems.
///
/// ### Params
///
/// * `mat_a` - Matrix A of the system
/// * `mat_b` - Matrix B of the system
/// * `mat_c` - Matrix C of the system
/// * `max_iter` - Maximum number of iterations
/// * `tolerance` - Tolerance parameter
///
/// Returns
///
/// The matrix X
fn sylvester_solver_iterative(
    mat_a: &MatRef<f64>,
    mat_b: &MatRef<f64>,
    mat_c: &MatRef<f64>,
    max_iter: usize,
    tolerance: f64,
) -> Mat<f64> {
    let m = mat_a.nrows();
    let n = mat_b.ncols();

    // Initial guess
    let mut x = mat_c.to_owned();
    let mut x_new = Mat::zeros(m, n);
    let mut residual = Mat::zeros(m, n);

    // Adaptive alpha parameters
    let mut alpha: f64 = 0.5;
    let alpha_min = 0.01;
    let alpha_max = 1.0;
    let mut prev_residual_norm = f64::INFINITY;

    let c_norm = mat_c.norm_l2();
    let rel_tolerance = tolerance * c_norm.max(1.0);

    for iter in 0..max_iter {
        // x_new = C - alpha * (A*x + x*B)
        let ax = mat_a * &x;
        let xb = &x * mat_b;

        residual.copy_from(&mat_c);
        residual -= &ax;
        residual -= &xb;

        let residual_norm = residual.norm_l2();

        // Check convergence with relative tolerance
        if residual_norm < rel_tolerance {
            break;
        }

        if iter > 0 {
            if residual_norm < prev_residual_norm {
                // Good progress, increase step size slightly
                alpha = (alpha * 1.1).min(alpha_max);
            } else {
                // Poor progress, reduce step size
                alpha = (alpha * 0.5).max(alpha_min);
            }
        }

        x_new.copy_from(&x);

        x_new.add_by_ref(&(residual.as_ref() * alpha));

        std::mem::swap(&mut x, &mut x_new);
        prev_residual_norm = residual_norm;

        // Early termination for very slow convergence
        if iter > 10 && residual_norm > 0.99 * prev_residual_norm {
            break;
        }
    }

    x
}

/// Direct version for small matrices
///
/// Uses partial LU decomposition to solve: `AX + XB = C`. Slow for large
/// matrix systems.
///
/// ### Params
///
/// * `mat_a` - Matrix A of the system
/// * `mat_b` - Matrix B of the system
/// * `mat_c` - Matrix C of the system
///
/// ### Returns
///
/// The matrix X
fn sylvester_solver_direct(
    mat_a: &MatRef<f64>,
    mat_b: &MatRef<f64>,
    mat_c: &MatRef<f64>,
) -> Mat<f64> {
    let m = mat_a.nrows();
    let n = mat_b.ncols();
    let mn = m * n;

    let mut coeff_matrix: Mat<f64> = Mat::zeros(mn, mn);

    // Build coefficient matrix
    for i in 0..m {
        for j in 0..n {
            let row_idx = i * n + j;

            // A part: (I ⊗ A)
            for k in 0..m {
                let col_idx = k * n + j;
                coeff_matrix[(row_idx, col_idx)] = mat_a[(i, k)];
            }

            // B^T part: (B^T ⊗ I)
            for l in 0..n {
                let col_idx = i * n + l;
                coeff_matrix[(row_idx, col_idx)] += mat_b[(l, j)];
            }
        }
    }

    // Vectorise C
    let mut c_vec: Mat<f64> = Mat::zeros(mn, 1);
    for i in 0..m {
        for j in 0..n {
            c_vec[(i * n + j, 0)] = mat_c[(i, j)];
        }
    }

    let lu = PartialPivLu::new(coeff_matrix.as_ref());
    let solved = lu.solve(&c_vec);

    // Reshape
    let mut res = Mat::zeros(m, n);
    for i in 0..m {
        for j in 0..n {
            res[(i, j)] = solved[(i * n + j, 0)];
        }
    }

    res
}

///////////
// Tests //
///////////

#[cfg(test)]
mod tests {
    use super::*;
    use faer::mat;

    #[test]
    fn test_sylvester_solver() {
        let a = mat![[1.0, 2.0], [0.0, 3.0]];
        let b = mat![[4.0, 1.0], [0.0, 2.0]];
        let c = mat![[1.0, 2.0], [3.0, 4.0]];

        let x = sylvester_solver(&a.as_ref(), &b.as_ref(), &c.as_ref());

        // Verify: AX + XB should equal C
        let result = &a * &x + &x * &b;

        // Check if close to C (allowing for numerical errors)
        for i in 0..c.nrows() {
            for j in 0..c.ncols() {
                let diff = result[(i, j)] - c[(i, j)].abs();
                assert!(diff < 1e-10, "Solution verification failed");
            }
        }
    }
}
