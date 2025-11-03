use arpack_ng::{eigenvectors, Which};
use num_complex::Complex64;
use std::ops::{Add, Mul};

use crate::core::data::sparse_structures::*;

/// Compute largest eigenvalues and eigenvectors using ARPACK-NG
///
/// ### Params
///
/// * `matrix` - Sparse matrix in CSR format
/// * `n_components` - Number of eigenpairs to compute
///
/// ### Returns
///
/// (eigenvalues, eigenvectors) where eigenvectors[i][j] is element j of
/// eigenvector i
pub fn compute_largest_eigenpairs_arpack<T>(
    matrix: &CompressedSparseData<T>,
    n_components: usize,
) -> (Vec<f32>, Vec<Vec<f32>>)
where
    T: Clone + Default + Into<f64> + Sync + Add + PartialEq + Mul,
{
    let n = matrix.shape.0;
    let csc = match matrix.cs_type {
        CompressedSparseFormat::Csc => matrix.clone(),
        CompressedSparseFormat::Csr => matrix.transform(),
    };

    eprintln!(
        "Matrix shape: {:?}, n_components: {}",
        matrix.shape, n_components
    );

    let ncv = (3 * n_components).min(n); // Lanczos vectors
    let maxiter = 300;

    let csc_f64: Vec<f64> = csc.data.iter().map(|v| v.clone().into()).collect();
    let csc_indices = csc.indices.clone();
    let csc_indptr = csc.indptr.clone();

    let (evals, evecs) = eigenvectors(
        move |x, mut y| {
            y.fill(Complex64::new(0.0, 0.0));
            for j in 0..n {
                if j + 1 >= csc_indptr.len() {
                    eprintln!("ERROR: j={} exceeds indptr length={}", j, csc_indptr.len());
                    continue;
                }
                for idx in csc_indptr[j]..csc_indptr[j + 1] {
                    if idx >= csc_indices.len() || idx >= csc_f64.len() {
                        eprintln!("ERROR: idx={} out of bounds", idx);
                        continue;
                    }
                    let i = csc_indices[idx];
                    if i >= n {
                        eprintln!("ERROR: i={} >= n={}", i, n);
                        continue;
                    }
                    let val = Complex64::new(csc_f64[idx], 0.0);
                    y[i] += val * x[j];
                }
            }
        },
        n,
        &Which::LargestMagnitude,
        n_components,
        ncv,
        maxiter,
    )
    .unwrap();

    let eigenvalues: Vec<f32> = evals.iter().map(|c| c.re as f32).collect();

    // Convert eigenvectors: evecs is shape (n, n_components)
    let mut eigenvectors = vec![vec![0.0f32; n_components]; n];
    for i in 0..n {
        for j in 0..n_components {
            eigenvectors[i][j] = evecs[[i, j]].re as f32;
        }
    }

    (eigenvalues, eigenvectors)
}
