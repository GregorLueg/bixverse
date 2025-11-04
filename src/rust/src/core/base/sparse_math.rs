use faer::traits::{ComplexField, RealField};
use faer::Mat;
use num_traits::Float;
use num_traits::ToPrimitive;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::iter::*;
use std::iter::Sum;
use std::ops::{Add, Mul};

use crate::assert_same_len;
use crate::core::data::sparse_structures::*;

/////////////////////////////////////
// Lanczos Eigenvalue calculations //
/////////////////////////////////////

/// Helper function for dot product of two vectors
///
/// ### Params
///
/// * `a` - Vector a
/// * `b` - Vector b
///
/// ### Returns
///
/// Dot product of the two vectors
fn dot<T>(a: &[T], b: &[T]) -> T
where
    T: Float + Send + Sync + Sum,
{
    assert_same_len!(a, b);
    a.par_iter().zip(b).map(|(&x, &y)| x * y).sum()
}

/// Helper function to normalise a vector
///
/// ### Params
///
/// * `v` - Initial vector
///
/// ### Returns
///
/// Normalised dot product of the vector `v`
fn norm<T>(v: &[T]) -> T
where
    T: Float + Send + Sync + Sum,
{
    dot(v, v).sqrt()
}

/// Helper function to normalise a vector
///
/// ### Params
///
/// * `v` - Mutable reference of the vector to normalise
fn normalise<T>(v: &mut [T])
where
    T: Float + Send + Sync + Sum,
{
    let n = norm(v);
    v.par_iter_mut().for_each(|x| *x = *x / n);
}

/// Helper function to calculate eigenvalues
///
/// ### Params
///
/// * `alpha` - alpha vector
/// * `beta` - beta vector
///
/// ### Returns
///
/// Tuple of `(eigenvectors, eigenvalues)`
fn tridiag_eig<T>(alpha: &[T], beta: &[T]) -> (Vec<T::Real>, Mat<T>)
where
    T: ComplexField + Copy + RealField,
{
    let n = alpha.len();
    let mut t = Mat::<T>::zeros(n, n);

    for i in 0..n {
        t[(i, i)] = alpha[i];
        if i < n - 1 {
            t[(i, i + 1)] = beta[i];
            t[(i + 1, i)] = beta[i];
        }
    }

    let eig = t.self_adjoint_eigen(faer::Side::Lower).unwrap();
    let evals = eig.S().column_vector().iter().copied().collect();
    let evecs = eig.U().to_owned();

    (evals, evecs)
}

/// Compute largest eigenvalues and eigenvectors using Lanczos
///
/// ### Params
///
/// * `matrix` - Sparse matrix in CSR format
/// * `n_components` - Number of eigenpairs to compute
/// * `seed` - For reproducibility
///
/// ### Returns
///
/// (eigenvalues, eigenvectors) where eigenvectors[i][j] is element j of
/// eigenvector i
pub fn compute_largest_eigenpairs_lanczos<T>(
    matrix: &CompressedSparseData<T>,
    n_components: usize,
    seed: u64,
) -> (Vec<f32>, Vec<Vec<f32>>)
where
    T: Clone + Default + Into<f64> + Sync + Add + PartialEq + Mul,
{
    let n = matrix.shape.0;
    let n_iter = (n_components * 2 + 10).max(n_components).min(n);

    // Convert to CSR for efficient row access
    let csr = match matrix.cs_type {
        CompressedSparseFormat::Csr => matrix.clone(),
        CompressedSparseFormat::Csc => matrix.transform(),
    };

    let data_f64: Vec<f64> = csr.data.iter().map(|v| v.clone().into()).collect();

    let matvec = |x: &[f64], y: &mut [f64]| {
        y.fill(0.0);
        for i in 0..n {
            for idx in csr.indptr[i]..csr.indptr[i + 1] {
                let j = csr.indices[idx];
                y[i] += data_f64[idx] * x[j];
            }
        }
    };

    // Lanczos iteration
    let mut v = vec![0.0; n];
    let mut v_old = vec![0.0; n];
    let mut w = vec![0.0; n];
    let mut v_matrix = vec![vec![0.0; n]; n_iter];

    let mut rng = StdRng::seed_from_u64(seed);

    for i in 0..n {
        v[i] = rng.random::<f64>() - 0.5;
    }
    normalise(&mut v);

    let mut alpha = vec![0.0; n_iter];
    let mut beta = vec![0.0; n_iter];

    for j in 0..n_iter {
        v_matrix[j].copy_from_slice(&v);

        matvec(&v, &mut w);
        alpha[j] = dot(&w, &v);

        // w = w - alpha[j]*v - beta[j-1]*v_old
        for i in 0..n {
            w[i] -= alpha[j] * v[i];
            if j > 0 {
                w[i] -= beta[j - 1] * v_old[i];
            }
        }

        beta[j] = norm(&w);
        if beta[j] < 1e-12 {
            break;
        }

        v_old.copy_from_slice(&v);
        v.copy_from_slice(&w);
        normalise(&mut v);
    }

    let (evals, evecs) = tridiag_eig(&alpha[..n_iter], &beta[..n_iter - 1]);

    let mut indices: Vec<usize> = (0..evals.len()).collect();
    indices.sort_by(|&i, &j| evals[j].partial_cmp(&evals[i]).unwrap());

    let mut largest_evals: Vec<f32> = Vec::with_capacity(n_components);
    let mut largest_evecs: Vec<Vec<f32>> = Vec::with_capacity(n_components);

    for &idx in indices.iter().take(n_components) {
        // Transform eigenvector back to original space: v_original = V * v_tridiag
        let mut evec = vec![0.0; n];
        for i in 0..n {
            for j in 0..n_iter {
                evec[i] += v_matrix[j][i] * evecs[(j, idx)].to_f64().unwrap();
            }
        }

        // Normalise the transformed eigenvector... Really should do this...
        let norm: f64 = evec.iter().map(|x| x * x).sum::<f64>().sqrt();
        for x in &mut evec {
            *x /= norm;
        }

        largest_evals.push(evals[idx].to_f64().unwrap() as f32);
        largest_evecs.push(evec.iter().map(|&x| x as f32).collect());
    }

    let mut transposed = vec![vec![0.0f32; n_components]; n];
    for comp_idx in 0..n_components {
        for point_idx in 0..n {
            transposed[point_idx][comp_idx] = largest_evecs[comp_idx][point_idx];
        }
    }

    (largest_evals, transposed)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lanczos_simple_matrix() {
        // Create symmetric matrix [[2, 1], [1, 2]]
        // Known eigenvalues: 3.0, 1.0
        let data = vec![2.0, 1.0, 1.0, 2.0];
        let indices = vec![0, 1, 0, 1];
        let indptr = vec![0, 2, 4];

        let matrix = CompressedSparseData {
            data: data.clone(),
            indices,
            indptr,
            shape: (2, 2),
            cs_type: CompressedSparseFormat::Csr,
            data_2: Some(data.to_vec()),
        };

        let (evals, evecs) = compute_largest_eigenpairs_lanczos(&matrix, 2, 42);

        println!("Eigenvalues: {:?}", evals);

        // Check largest eigenvalue is approximately 3.0
        assert!(
            (evals[0] - 3.0).abs() < 0.01,
            "Expected largest eigenvalue ~3.0, got {}",
            evals[0]
        );
        assert!(
            (evals[1] - 1.0).abs() < 0.01,
            "Expected second eigenvalue ~1.0, got {}",
            evals[1]
        );

        // Check eigenvector dimensions
        assert_eq!(evecs.len(), 2, "Should have 2 eigenvectors");
        assert_eq!(evecs[0].len(), 2, "Each eigenvector should have 2 elements");
    }
}
