use extendr_api::prelude::*;
use faer::Mat;
use rayon::iter::*;
use rustc_hash::FxHashSet;

use crate::core::base::cors_similarity::*;
use crate::utils::general::{string_vec_to_set, upper_triangle_indices};
use crate::utils::r_rust_interface::*;

/// Calculate the column-wise co-variance.
///
/// @description Calculates the co-variance of the columns.
/// WARNING! Incorrect use can cause kernel crashes. Wrapper around the Rust
/// functions with type checks are provided in the package.
///
/// @param x R matrix with doubles.
///
/// @returns The co-variance matrix.
///
/// @export
#[extendr]
fn rs_covariance(x: RMatrix<f64>) -> extendr_api::RArray<f64, [usize; 2]> {
    let mat = r_matrix_to_faer(&x);
    let covar = column_cov(&mat);

    faer_to_r_matrix(covar.as_ref())
}

/// Calculate the column wise correlations.
///
/// @description Calculates the correlation matrix of the columns.
/// WARNING! Incorrect use can cause kernel crashes. Wrapper around the Rust
/// functions with type checks are provided in the package.
///
/// @param x R matrix with doubles.
/// @param spearman Shall the Spearman correlation be calculated instead of
/// Pearson.
///
/// @returns The correlation matrix.
///
/// @export
#[extendr]
fn rs_cor(x: RMatrix<f64>, spearman: bool) -> extendr_api::RArray<f64, [usize; 2]> {
    let mat = r_matrix_to_faer(&x);

    let cor = column_cor(&mat, spearman);

    faer_to_r_matrix(cor.as_ref())
}

/// Calculate the column wise cosine similarities
///
/// @description Calculates the cosyne similarity matrix of the columns.
///
/// @param x R matrix with doubles.
///
/// @returns The correlation matrix.
///
/// @export
#[extendr]
fn rs_cos(x: RMatrix<f64>) -> extendr_api::RArray<f64, [usize; 2]> {
    let mat = r_matrix_to_faer(&x);

    let cos = column_cos(&mat);

    faer_to_r_matrix(cos.as_ref())
}

/// Calculate the column wise correlations.
///
/// @description Calculates the correlation between the columns of two matrices.
/// The number of rows need to be the same!
///
/// @param x R matrix with doubles.
/// @param y R matrix with doubles.
/// @param spearman Shall the Spearman correlation be calculated instead of
/// Pearson.
///
/// @returns The correlation matrix.
///
/// @export
#[extendr]
fn rs_cor2(
    x: RMatrix<f64>,
    y: RMatrix<f64>,
    spearman: bool,
) -> extendr_api::RArray<f64, [usize; 2]> {
    let x = r_matrix_to_faer(&x);
    let y = r_matrix_to_faer(&y);

    let cor = cor(&x, &y, spearman);

    faer_to_r_matrix(cor.as_ref())
}

/// Calculates the correlation matrix from the co-variance matrix
///
/// @description Calculates the correlation matrix from a co-variance
/// matrix
///
/// @param x R matrix with doubles that is the co-variance matrix
///
/// @returns The correlation matrix.
///
/// @export
#[extendr]
fn rs_cov2cor(x: RMatrix<f64>) -> extendr_api::RArray<f64, [usize; 2]> {
    let mat = r_matrix_to_faer(&x);

    let cor = cov2cor(mat);

    faer_to_r_matrix(cor.as_ref())
}

/// Calculates the mutual information matrix
///
/// @description Calculates the mutual information across all columns in the
/// data.
///
/// @param x R matrix with doubles for which to calculate the mutual information
/// @param n_bins Optional integer. Number of bins to use. If `NULL` is provided
/// the function will default to `sqrt(nrows(x))`.
/// @param strategy String. Binning strategy One of
/// `c("equal_width", "equal_freq")`.
/// @param normalise Boolean. Shall the normalised mutual information be
/// calculated via joint entropy.
///
/// @returns The mutual information matrix.
///
/// @export
#[extendr]
fn rs_mutual_info(
    x: RMatrix<f64>,
    n_bins: Option<usize>,
    strategy: String,
    normalise: bool,
) -> extendr_api::Result<extendr_api::RArray<f64, [usize; 2]>> {
    let mat = r_matrix_to_faer(&x);

    let mi_mat = column_mutual_information(&mat, n_bins, normalise, &strategy)?;

    Ok(faer_to_r_matrix(mi_mat.as_ref()))
}

/// Calculates the point wise mutual information
///
/// @description Calculates the pointwise mutual information (can be also
/// normalised) across all columns of the data. This can be used to identify
/// (dis)similar samples based on Boolean characteristics.
///
/// @param x Logical matrix. The columns represent features and the rows
/// represent samples
/// @param normalise Shall the normalised pointwise mutual information be
/// returned.
///
/// @returns The (normalised) pointwise mutual information matrix.
///
/// @export
#[extendr]
fn rs_pointwise_mutual_info(
    x: RMatrix<Rbool>,
    normalise: bool,
) -> extendr_api::RArray<f64, [usize; 2]> {
    let data = r_matrix_to_vec_bool(&x);

    let npmi_mat = calc_pmi(&data, normalise);

    faer_to_r_matrix(npmi_mat.as_ref())
}

/// Calculate the pairwise column distance in a matrix
///
/// @description
/// This function allows to calculate pairwise between all columns the specified
/// distance metric.
///
/// @param x Numerical matrix. The matrix for which to calculate the pairwise
/// column distances.
/// @param distance_type String. One of
/// `c("euclidean", "manhattan", "canberra", "cosine")`.
///
/// @return The calculated distance matrix
///
/// @export
#[extendr]
fn rs_dist(
    x: RMatrix<f64>,
    distance_type: String,
) -> extendr_api::Result<extendr_api::RArray<f64, [usize; 2]>> {
    let data = &r_matrix_to_faer(&x);

    let dist_type = parse_distance_type(&distance_type)
        .ok_or_else(|| format!("Invalid Distance type: {}", distance_type))?;

    let res = match dist_type {
        DistanceType::L2Norm => column_pairwise_l2_norm(data),
        DistanceType::L1Norm => column_pairwise_l1_norm(data),
        DistanceType::Cosine => column_pairwise_cosine_dist(data),
        DistanceType::Canberra => column_pairwise_canberra_dist(data),
    };

    Ok(faer_to_r_matrix(res.as_ref()))
}

/// Calculate the column wise correlations.
///
/// @description Calculates the correlation matrix of the columns. This function
/// will return the upper triangle. WARNING! Incorrect use can cause kernel
/// crashes. Wrapper around the Rust functions with type checks are provided in
/// the package.
///
/// @param x R matrix with doubles.
/// @param spearman Shall the Spearman correlation be calculated instead of
/// Pearson.
/// @param shift Shall a shift be applied to the matrix. 0 = the diagonal will
/// be included. 1 = the diagonal will not be included.
///
/// @returns The upper triangle of the correlation matrix iterating through the
/// rows, shifted by one (the diagonal will not be returned).
///
/// @export
#[extendr]
fn rs_cor_upper_triangle(x: RMatrix<f64>, spearman: bool, shift: usize) -> Vec<f64> {
    // Calculate the correlations
    let mat = r_matrix_to_faer(&x);
    let cor = column_cor(&mat, spearman);
    let upper_triangle_indices = upper_triangle_indices(mat.ncols(), shift);
    let mut cor_flat = Vec::new();
    for (&r, &c) in upper_triangle_indices
        .0
        .iter()
        .zip(upper_triangle_indices.1.iter())
    {
        cor_flat.push(*cor.get(r, c));
    }

    cor_flat
}

/// Set similarities
///
/// This function calculates the Jaccard or similarity index between a two given
/// string vector and a  of other string vectors.
///
/// @param s_1 The String vector against which to calculate the set similarities.
/// @param s_2 The String vector against which to calculate the set similarities.
/// @param overlap_coefficient Boolean. Use the overlap coefficient instead of the Jaccard similarity be calculated.
///
/// @export
#[extendr]
fn rs_set_similarity(s_1: Vec<String>, s_2: Vec<String>, overlap_coefficient: bool) -> f64 {
    let mut s_hash1 = FxHashSet::default();
    let mut s_hash2 = FxHashSet::default();

    for item in &s_1 {
        s_hash1.insert(item);
    }
    for item in &s_2 {
        s_hash2.insert(item);
    }

    set_similarity(&s_hash1, &s_hash2, overlap_coefficient)
}

/// Set similarities over two list
///
/// @description
/// This function calculates the Jaccard or similarity index between two lists.
///
/// @param s_1_list R list. The first list of string elements you want to
/// compare against.
/// @param s_2_list R list. The second list of string elements you want to
/// compare against.
/// @param overlap_coefficient Boolean. Use the overlap coefficient instead of
/// the Jaccard similarity be calculated.
///
/// @return A matrix of the Jaccard similarities between the elements. The rows
/// represent `s_1_list` and the column `s_2_list`.
///
/// @export
#[extendr]
fn rs_set_similarity_list2(
    s_1_list: List,
    s_2_list: List,
    overlap_coefficient: bool,
) -> extendr_api::Result<RArray<f64, [usize; 2]>> {
    let s1_vec = r_list_to_str_vec(s_1_list)?;
    let s2_vec = r_list_to_str_vec(s_2_list)?;

    let s_hash1: Vec<FxHashSet<&String>> = s1_vec.iter().map(|s| string_vec_to_set(s)).collect();
    let s_hash2: Vec<FxHashSet<&String>> = s2_vec.iter().map(|s| string_vec_to_set(s)).collect();

    let data: Vec<f64> = s_hash1
        .par_iter()
        .flat_map(|s1| {
            s_hash2
                .iter()
                .map(|s2| set_similarity(s1, s2, overlap_coefficient))
                .collect::<Vec<_>>()
        })
        .collect();

    let matrix = Mat::from_fn(s_hash1.len(), s_hash2.len(), |i, j| {
        data[i * s_hash2.len() + j]
    });

    Ok(faer_to_r_matrix(matrix.as_ref()))
}

/// Set similarities over one list
///
/// @description This function calculates the set similarity via Jaccard or
/// overlap coefficient across all permutations of one list.
///
/// @param list A named R list.
/// @param overlap_coefficient Boolean. Use the overlap coefficient instead of
/// the Jaccard similarity be calculated.
///
/// @return A list with the following items:
/// \itemize{
///     \item from - Name of element i
///     \item to - Name of element j
///     \item sim - Similarity between the two elements
/// }
///
/// @export
#[extendr]
fn rs_set_similarity_list(list: List, overlap_coefficient: bool) -> extendr_api::Result<List> {
    let btree = r_list_to_btree_set(list)?;
    let names: Vec<&String> = btree.keys().collect();

    let mut name_i = Vec::new();
    let mut name_j = Vec::new();
    let mut sim = Vec::new();

    for i in 0..names.len() {
        for j in (i + 1)..names.len() {
            name_i.push(names[i].clone());
            name_j.push(names[j].clone());
            let val_i = btree.get(names[i]).unwrap().iter().collect();
            let val_j = btree.get(names[j]).unwrap().iter().collect();
            sim.push(set_similarity(&val_i, &val_j, overlap_coefficient));
        }
    }

    Ok(list!(from = name_i, to = name_j, sim = sim))
}

extendr_module! {
  mod r_cors_similarity;
  fn rs_covariance;
  fn rs_cor;
  fn rs_cos;
  fn rs_cor2;
  fn rs_cov2cor;
  fn rs_dist;
  fn rs_mutual_info;
  fn rs_pointwise_mutual_info;
  fn rs_cor_upper_triangle;
  fn rs_set_similarity;
  fn rs_set_similarity_list;
  fn rs_set_similarity_list2;
}
