use faer::{Mat, MatRef};
use rayon::iter::*;

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
