use faer::{
    linalg::solvers::{PartialPivLu, Solve},
    traits::AddByRef,
    Mat, MatRef,
};
use rand::prelude::*;
use rand_distr::Normal;
use rayon::iter::*;

use crate::utils::general::*;
use crate::utils::utils_stats::*;

use crate::{assert_nrows, assert_same_dims, assert_symmetric_mat};

//////////////////////////////
// ENUMS, TYPES, STRUCTURES //
//////////////////////////////

/// Structure for random SVD results
///
/// ### Fields
///
/// * `u` - Matrix u of the SVD decomposition
/// * `v` - Matrix v of the SVD decomposition
/// * `s` - Eigen vectors of the SVD decomposition
#[derive(Clone, Debug)]
pub struct RandomSvdResults {
    pub u: faer::Mat<f64>,
    pub v: faer::Mat<f64>,
    pub s: Vec<f64>,
}

/// Structure for DiffCor results
///
/// ### Fields
///
/// * `r_a` - Correlation coefficients of a
/// * `r_b` - Correlation coefficients of b
/// * `z_score` - Z-scores of the differential correlation
/// * `p_vals` - Calculated p-values from the Z-scores
#[derive(Clone, Debug)]
pub struct DiffCorRes {
    pub r_a: Vec<f64>,
    pub r_b: Vec<f64>,
    pub z_score: Vec<f64>,
    pub p_vals: Vec<f64>,
}

//////////////////////////////
// SCALING, COVAR, COR, PCA //
//////////////////////////////

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
/// * `mat` - The matrix for which to calculate the column-wise standard
///           deviations
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
fn normalise_matrix_col_l2(mat: &MatRef<f64>) -> Mat<f64> {
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

/// Calculate the co-variance
///
/// ### Params
///
/// * `mat` - The matrix for which to calculate the co-variance. Assumes that
///           features are columns.
///
/// ### Returns
///
/// The resulting co-variance matrix.
pub fn column_covariance(mat: &MatRef<f64>) -> Mat<f64> {
    assert_symmetric_mat!(mat);

    let n_rows = mat.nrows();
    let centered = scale_matrix_col(mat, false);
    let covariance = (centered.transpose() * &centered) / (n_rows - 1) as f64;

    covariance
}

/// Calculate the cosine similarity
///
/// ### Params
///
/// * `mat` - The matrix for which to calculate the cosine similarity. Assumes
///           that features are columns.
///
/// ### Returns
///
/// The resulting cosine similarity matrix
pub fn column_cosine(mat: &MatRef<f64>) -> Mat<f64> {
    let normalised = normalise_matrix_col_l2(mat);

    normalised.transpose() * &normalised
}

/// Calculate the correlation matrix
///
/// ### Params
///
/// * `mat` - The matrix for which to calculate the correlation matrix. Assumes
///           that features are columns.
/// * `spearman` - Shall Spearman correlation be used.
///
/// ### Returns
///
/// The resulting correlation matrix.
pub fn column_correlation(mat: &MatRef<f64>, spearman: bool) -> Mat<f64> {
    let mat = if spearman {
        rank_matrix_col(mat)
    } else {
        mat.to_owned()
    };

    let scaled = scale_matrix_col(&mat.as_ref(), true);

    let nrow = scaled.nrows() as f64;

    let cor = scaled.transpose() * &scaled / (nrow - 1_f64);

    cor
}

/// Calculates the correlation between two matrices
///
/// The two matrices need to have the same number of rows, otherwise the function
/// panics
///
/// ### Params
///
/// * `mat_a` - The first matrix.
/// * `mat_b` - The second matrix.
/// * `spearman` - Shall Spearman correlation be used.
///
/// ### Returns
///
/// The resulting correlation between the samples of the two matrices
pub fn cor(mat_a: &MatRef<f64>, mat_b: &MatRef<f64>, spearman: bool) -> Mat<f64> {
    assert_nrows!(mat_a, mat_b);

    let nrow = mat_a.nrows() as f64;

    let mat_a = if spearman {
        rank_matrix_col(mat_a)
    } else {
        mat_a.to_owned()
    };

    let mat_b = if spearman {
        rank_matrix_col(mat_b)
    } else {
        mat_b.to_owned()
    };

    let mat_a = scale_matrix_col(&mat_a.as_ref(), true);
    let mat_b = scale_matrix_col(&mat_b.as_ref(), true);

    let cor = mat_a.transpose() * &mat_b / (nrow - 1_f64);

    cor
}

/// Calculate the correlation matrix from the co-variance matrix
///
/// ### Params
///
/// * `mat` - The co-variance matrix
///
/// ### Returns
///
/// The resulting correlation matrix.
pub fn cov2cor(mat: MatRef<f64>) -> Mat<f64> {
    assert_symmetric_mat!(mat);

    let n = mat.nrows();
    let mut result = mat.to_owned();

    let inv_sqrt_diag: Vec<f64> = (0..n).map(|i| 1.0 / mat.get(i, i).sqrt()).collect();

    for i in 0..n {
        for j in 0..n {
            result[(i, j)] = mat.get(i, j) * inv_sqrt_diag[i] * inv_sqrt_diag[j];
        }
    }

    result
}

/// Get the eigenvalues and vectors from a covar or cor matrix
///
/// Function will panic if the matrix is not symmetric
///
/// ### Params
///
/// * `matrix` - The correlation or co-variance matrix
/// * `top_n` - How many of the top eigen vectors and values to return.
///
/// ### Returns
///
/// A vector of tuples corresponding to the top eigen pairs.
pub fn get_top_eigenvalues(matrix: &Mat<f64>, top_n: usize) -> Vec<(f64, Vec<f64>)> {
    // Ensure the matrix is square
    assert_symmetric_mat!(matrix);

    let eigendecomp = matrix.eigen().unwrap();

    let s = eigendecomp.S();
    let u = eigendecomp.U();

    // Extract the real part of the eigenvalues and vectors
    let mut eigenpairs = s
        .column_vector()
        .iter()
        .zip(u.col_iter())
        .map(|(l, v)| {
            let l_real = l.re;
            let v_real = v.iter().map(|v_i| v_i.re).collect::<Vec<f64>>();
            (l_real, v_real)
        })
        .collect::<Vec<(f64, Vec<f64>)>>();

    // Sort and return Top N
    eigenpairs.sort_by(|a, b| b.0.total_cmp(&a.0));

    let res: Vec<(f64, Vec<f64>)> = eigenpairs.into_iter().take(top_n).collect();

    res
}

/// Randomised SVD
///
/// ### Params
///
/// * `x` - The matrix on which to apply the randomised SVD.
/// * `rank` - The target rank of the approximation (number of singular values,
///            vectors to compute).
/// * `seed` - Random seed for reproducible results.
/// * `oversampling` - Additional samples beyond the target rank to improve accuracy.
///                    Defaults to 10 if not specified.
/// * `n_power_iter` - Number of power iterations to perform for better approximation quality.
///                    More iterations generally improve accuracy but increase computation time.
///                    Defaults to 2 if not specified.
///
/// ### Returns
///
/// The randomised SVD results in form of `RandomSvdResults`.
///
/// ### Algorithm Details
///
/// 1. Generate a random Gaussian matrix Ω of size n × (rank + oversampling)
/// 2. Compute Y = X * Ω to capture the range of X
/// 3. Orthogonalize Y using QR decomposition to get Q
/// 4. Apply power iterations: for each iteration, compute Z = X^T * Q, then Q = QR(X * Z)
/// 5. Form B = Q^T * X and compute its SVD
/// 6. Reconstruct the final SVD: U = Q * U_B, V = V_B, S = S_B
pub fn randomised_svd(
    x: MatRef<f64>,
    rank: usize,
    seed: usize,
    oversampling: Option<usize>,
    n_power_iter: Option<usize>,
) -> RandomSvdResults {
    let ncol = x.ncols();
    let nrow = x.nrows();

    // Oversampling for better accuracy
    let os = oversampling.unwrap_or(10);
    let sample_size = (rank + os).min(ncol.min(nrow));
    let n_iter = n_power_iter.unwrap_or(2);

    // Create a random matrix
    let mut rng = StdRng::seed_from_u64(seed as u64);
    let normal = Normal::new(0.0, 1.0).unwrap();
    let omega = Mat::from_fn(ncol, sample_size, |_, _| normal.sample(&mut rng));

    // Multiply random matrix with original and use QR composition to get
    // low rank approximation of x
    let y = x * omega;

    let mut q = y.qr().compute_thin_Q();
    for _ in 0..n_iter {
        let z = x.transpose() * q;
        q = (x * z).qr().compute_thin_Q();
    }

    // Perform the SVD on the low-rank approximation
    let b = q.transpose() * x;
    let svd = b.thin_svd().unwrap();

    RandomSvdResults {
        u: q * svd.U(),
        v: svd.V().cloned(), // Use clone instead of manual copying
        s: svd.S().column_vector().iter().copied().collect(),
    }
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
// Other //
///////////

/// Calculate differential correlations
///
/// The function will panic if the two correlation matrices are not symmetric
/// and do not have the same dimensions.
///
/// ### Params
///
/// * `mat_a` - The first correlation matrix.
/// * `mat_b` - The second correlation matrix.
/// * `no_sample_a` - Number of samples that were present to calculate mat_a.
/// * `no_sample_b` - Number of samples that were present to calculate mat_b.
/// * `spearman` - Was Spearman correlation used.
///
/// ### Returns
///
/// The resulting differential correlation results as a structure.
pub fn calculate_diff_correlation(
    mat_a: &Mat<f64>,
    mat_b: &Mat<f64>,
    no_sample_a: usize,
    no_sample_b: usize,
    spearman: bool,
) -> DiffCorRes {
    assert_symmetric_mat!(mat_a);
    assert_symmetric_mat!(mat_b);
    assert_same_dims!(mat_a, mat_b);

    let mut cors_a: Vec<f64> = Vec::new();
    let mut cors_b: Vec<f64> = Vec::new();

    let upper_triangle_indices = upper_triangle_indices(mat_a.ncols(), 1);

    for (&r, &c) in upper_triangle_indices
        .0
        .iter()
        .zip(upper_triangle_indices.1.iter())
    {
        cors_a.push(*mat_a.get(r, c));
        cors_b.push(*mat_b.get(r, c));
    }

    // Maybe save the original correlations... Note to myself.
    let original_cor_a = cors_a.to_vec();
    let original_cor_b = cors_b.to_vec();

    cors_a.par_iter_mut().for_each(|x| *x = x.atanh());
    cors_b.par_iter_mut().for_each(|x| *x = x.atanh());

    // Constant will depend on if Spearman or Pearson
    let constant = if spearman { 1.06 } else { 1.0 };
    let denominator =
        ((constant / (no_sample_a as f64 - 3.0)) + (constant / (no_sample_b as f64 - 3.0))).sqrt();

    let z_scores: Vec<f64> = cors_a
        .par_iter()
        .zip(cors_b.par_iter())
        .map(|(a, b)| (a - b) / denominator)
        .collect();

    let p_values = z_scores_to_pval(&z_scores);

    DiffCorRes {
        r_a: original_cor_a,
        r_b: original_cor_b,
        z_score: z_scores,
        p_vals: p_values,
    }
}

/// Test different epsilons over a distance vector
///
/// ### Params
///
/// * `dist` - The distance vector on which to apply the specified RBF function.
///            Assumes that these are the values of upper triangle of the distance
///            matrix.
/// * `epsilons` - Vector of epsilons to test.
/// * `n` - Original dimensions of the distance matrix from which `dist` was
///         derived.
/// * `shift` - Was a shift applied during the generation of the vector, i.e., was
///             the diagonal included or not.
/// * `rbf_type` - Which RBF function to apply on the distance vector.
///
/// ### Returns
///
/// The column sums of the resulting adjacency matrices after application of the
/// RBF function to for example check if these are following power law distributions.
pub fn rbf_iterate_epsilons(
    dist: &[f64],
    epsilons: &[f64],
    n: usize,
    shift: usize,
    rbf_type: &str,
) -> Result<faer::Mat<f64>, String> {
    // Now specifying String as the error type
    let rbf_fun =
        parse_rbf_types(rbf_type).ok_or_else(|| format!("Invalid RBF function: {}", rbf_type))?;

    let k_res: Vec<Vec<f64>> = epsilons
        .par_iter()
        .map(|epsilon| {
            let affinity_adj = match rbf_fun {
                RbfType::Gaussian => rbf_gaussian(dist, epsilon),
                RbfType::Bump => rbf_bump(dist, epsilon),
                RbfType::InverseQuadratic => rbf_inverse_quadratic(dist, epsilon),
            };
            let affinity_adj_mat = upper_triangle_to_sym_faer(&affinity_adj, shift, n);
            col_sums(affinity_adj_mat.as_ref())
        })
        .collect();

    Ok(nested_vector_to_faer_mat(k_res, true))
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
