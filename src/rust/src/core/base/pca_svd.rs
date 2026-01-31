use faer::{Mat, MatRef};
use rand::prelude::*;
use rand_distr::Normal;

////////////////
// Structures //
////////////////

/// Structure for random SVD results
///
/// ### Fields
///
/// * `u` - Matrix u of the SVD decomposition
/// * `v` - Matrix v of the SVD decomposition
/// * `s` - Eigen vectors of the SVD decomposition
#[derive(Clone, Debug)]
pub struct RandomSvdResults<T> {
    pub u: faer::Mat<T>,
    pub v: faer::Mat<T>,
    pub s: Vec<T>,
}

///////////////
// Functions //
///////////////

/// Randomised SVD
///
/// ### Params
///
/// * `x` - The matrix on which to apply the randomised SVD.
/// * `rank` - The target rank of the approximation (number of singular values,
///   vectors to compute).
/// * `seed` - Random seed for reproducible results.
/// * `oversampling` - Additional samples beyond the target rank to improve
///   accuracy. Defaults to 10 if not specified.
/// * `n_power_iter` - Number of power iterations to perform for better
///   approximation quality. More iterations generally improve accuracy but
///   increase computation time. Defaults to 2 if not specified.
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
) -> RandomSvdResults<f64> {
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

/// Randomised SVD (f32)
///
/// The trait bound errors were crazy... Easy to just re-implement for f32.
///
/// ### Params
///
/// * `x` - The matrix on which to apply the randomised SVD.
/// * `rank` - The target rank of the approximation (number of singular values,
///   vectors to compute).
/// * `seed` - Random seed for reproducible results.
/// * `oversampling` - Additional samples beyond the target rank to improve accuracy.
///   Defaults to 10 if not specified.
/// * `n_power_iter` - Number of power iterations to perform for better approximation quality.
///   More iterations generally improve accuracy but increase computation time.
///   Defaults to 2 if not specified.
///
/// ### Returns
///
/// The randomised SVD results in form of `RandomSvdResults`.
pub fn randomised_svd_f32(
    x: MatRef<f32>,
    rank: usize,
    seed: usize,
    oversampling: Option<usize>,
    n_power_iter: Option<usize>,
) -> RandomSvdResults<f32> {
    let ncol = x.ncols();
    let nrow = x.nrows();

    // Oversampling for better accuracy
    let os = oversampling.unwrap_or(10);
    let sample_size = (rank + os).min(ncol.min(nrow));
    let n_iter = n_power_iter.unwrap_or(2);

    // Create a random matrix
    let mut rng = StdRng::seed_from_u64(seed as u64);
    let normal = Normal::new(0.0, 1.0).unwrap();
    let omega = Mat::from_fn(ncol, sample_size, |_, _| normal.sample(&mut rng) as f32);

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
