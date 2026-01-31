use faer::{Mat, MatRef};
use rustc_hash::FxHashSet;
use std::borrow::Borrow;
use std::hash::Hash;

use crate::assert_nrows;
use crate::core::base::utils::*;

///////////////////////
// Column matrix ops //
///////////////////////

/// Calculate the correlation matrix
///
/// ### Params
///
/// * `mat` - The matrix for which to calculate the correlation matrix. Assumes
///   that features are columns.
/// * `spearman` - Shall Spearman correlation be used.
///
/// ### Returns
///
/// The resulting correlation matrix.
pub fn column_cor(mat: &MatRef<f64>, spearman: bool) -> Mat<f64> {
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

//////////////////////
// Set similarities //
//////////////////////

/// Calculate the set similarity.
///
/// ### Params
///
/// * `s_1` - The first HashSet.
/// * `s_2` - The second HashSet.
/// * `overlap_coefficient` - Shall the overlap coefficient be returned or the
///   Jaccard similarity
///
/// ### Return
///
/// The Jaccard similarity or overlap coefficient.
pub fn set_similarity<T>(s_1: &FxHashSet<T>, s_2: &FxHashSet<T>, overlap_coefficient: bool) -> f64
where
    T: Borrow<String> + Hash + Eq,
{
    let i = s_1.intersection(s_2).count() as u64;
    let u = if overlap_coefficient {
        std::cmp::min(s_1.len(), s_2.len()) as u64
    } else {
        s_1.union(s_2).count() as u64
    };
    i as f64 / u as f64
}
