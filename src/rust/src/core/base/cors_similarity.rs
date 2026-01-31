use faer::{Mat, MatRef};
use rayon::iter::*;
use rustc_hash::FxHashSet;
use std::borrow::Borrow;
use std::hash::Hash;

use crate::assert_nrows;
use crate::core::base::utils::*;

//////////////////////////////
// ENUMS, TYPES, STRUCTURES //
//////////////////////////////

/// Distance type enum
#[derive(Debug, Clone)]
pub enum DistanceType {
    /// L1 norm, i.e., Manhattan distance
    L1Norm,
    /// L2 norm, i.e., Euclidean distance
    L2Norm,
    /// Cosine distance
    Cosine,
    /// Canberra distance
    Canberra,
}

////////////
// Params //
////////////

/// Parsing the distance type
///
/// ### Params
///
/// * `s` - string defining the distance type
///
/// ### Returns
///
/// The `DistanceType`.
pub fn parse_distance_type(s: &str) -> Option<DistanceType> {
    match s.to_lowercase().as_str() {
        "euclidean" => Some(DistanceType::L2Norm),
        "manhattan" => Some(DistanceType::L1Norm),
        "canberra" => Some(DistanceType::Canberra),
        "cosine" => Some(DistanceType::Cosine),
        _ => None,
    }
}

///////////////////////
// Column matrix ops //
///////////////////////

/// Calculate the cosine similarity between columns
///
/// ### Params
///
/// * `mat` - The matrix for which to calculate the cosine similarity.
///
/// ### Returns
///
/// The resulting cosine similarity matrix
pub fn column_cos(mat: &MatRef<f64>) -> Mat<f64> {
    let normalised = normalise_matrix_col_l2(mat);

    normalised.transpose() * &normalised
}

/// Calculate the cosine distance between columns
///
/// ### Params
///
/// * `mat` - The matrix for which to calculate the cosine distance.
///
/// ### Returns
///
/// The resulting cosine distance matrix
pub fn column_pairwise_cosine_dist(mat: &MatRef<f64>) -> Mat<f64> {
    let cosine_sim = column_cos(mat);
    let ncols = cosine_sim.ncols();
    let mut res: Mat<f64> = Mat::zeros(ncols, ncols);
    for i in 0..ncols {
        for j in 0..ncols {
            if i != j {
                res[(i, j)] = 1_f64 - cosine_sim.get(i, j).abs();
            }
        }
    }
    res
}

/// Calculate L2 Norm (Euclidean distance) between columns
///
/// ### Params
///
/// * `mat` - The matrix for which to calculate the L2 norm (Euclidean distance)
///   pairwise between all columns.
///
/// ### Returns
///
/// The distance matrix based on the L2 Norm.
pub fn column_pairwise_l2_norm(mat: &MatRef<f64>) -> Mat<f64> {
    let ncols = mat.ncols();

    let gram = mat.transpose() * mat;
    let mut col_norms_square = vec![0_f64; ncols];

    for i in 0..ncols {
        col_norms_square[i] = *gram.get(i, i);
    }

    let mut res: Mat<f64> = Mat::zeros(ncols, ncols);

    for i in 0..ncols {
        for j in 0..ncols {
            if i == j {
                res[(i, j)] = 0_f64;
            } else {
                let dist_sq = col_norms_square[i] + col_norms_square[j] - 2.0 * gram.get(i, j);
                let dist = dist_sq.max(0.0).sqrt();
                res[(i, j)] = dist;
            }
        }
    }

    res
}

/// Calculate L1 Norm (Manhatten distance) between columns
///
/// ### Params
///
/// * `mat` - The matrix for which to calculate the L1 norm (Manhatten distance)
///   pairwise between all columns.
///
/// ### Returns
///
/// The distance matrix based on the L1 Norm.
pub fn column_pairwise_l1_norm(mat: &MatRef<f64>) -> Mat<f64> {
    let (nrows, ncols) = mat.shape();
    let mut res: Mat<f64> = Mat::zeros(ncols, ncols);

    let pairs: Vec<(usize, usize)> = (0..ncols)
        .flat_map(|i| ((i + 1)..ncols).map(move |j| (i, j)))
        .collect();

    let results: Vec<(usize, usize, f64)> = pairs
        .par_iter()
        .map(|&(i, j)| {
            let col_i = mat.col(i);
            let col_j = mat.col(j);

            let mut sum_abs = 0_f64;

            // smid friendly code with unsafe
            let mut k = 0_usize;

            // 4 elements at once
            while k + 3 < nrows {
                unsafe {
                    sum_abs += (col_i.get_unchecked(k) - col_j.get_unchecked(k)).abs();
                    sum_abs += (col_i.get_unchecked(k + 1) - col_j.get_unchecked(k + 1)).abs();
                    sum_abs += (col_i.get_unchecked(k + 2) - col_j.get_unchecked(k + 2)).abs();
                    sum_abs += (col_i.get_unchecked(k + 3) - col_j.get_unchecked(k + 3)).abs();
                }
                k += 4;
            }

            // remaining elements
            while k < nrows {
                unsafe {
                    sum_abs += (col_i.get_unchecked(k) - col_j.get_unchecked(k)).abs();
                }
                k += 1;
            }

            (i, j, sum_abs)
        })
        .collect();

    for (i, j, dist) in results {
        res[(i, j)] = dist;
        res[(j, i)] = dist;
    }

    res
}

/// Calculate Canberra distance between columns
///
/// ### Params
///
/// * `mat` - The matrix for which to calculate the Canberra distance pairwise
///   between all columns.
///
/// ### Returns
///
/// The Canberra distance matrix between all columns.
pub fn column_pairwise_canberra_dist(mat: &MatRef<f64>) -> Mat<f64> {
    let (nrows, ncols) = mat.shape();
    let mut res: Mat<f64> = Mat::zeros(ncols, ncols);

    let pairs: Vec<(usize, usize)> = (0..ncols)
        .flat_map(|i| ((i + 1)..ncols).map(move |j| (i, j)))
        .collect();

    let results: Vec<(usize, usize, f64)> = pairs
        .par_iter()
        .map(|&(i, j)| {
            let col_i = mat.col(i);
            let col_j = mat.col(j);

            let mut sum_canberra = 0_f64;
            for k in 0..nrows {
                let val_i = col_i.get(k);
                let val_j = col_j.get(k);

                let abs_i = val_i.abs();
                let abs_j = val_j.abs();
                let denom = abs_i + abs_j;

                if denom > f64::EPSILON {
                    sum_canberra += (val_i - val_j).abs() / denom;
                }
            }

            (i, j, sum_canberra)
        })
        .collect();

    for (i, j, dist) in results {
        res[(i, j)] = dist;
        res[(j, i)] = dist;
    }

    res
}

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

////////////////////////////////////////
// Binary and other distance measures //
////////////////////////////////////////

/// Calculate Hamming distance between columns
///
/// ### Params
///
/// * `mat` - Integer matrix where categorical values are encoded as integers
///
/// ### Returns
///
/// The Hamming distance matrix with values in [0, 1]
pub fn column_pairwise_hamming_cat(mat: &MatRef<i32>) -> Mat<f64> {
    let (nrows, ncols) = mat.shape();
    let mut res = Mat::zeros(ncols, ncols);

    let pairs: Vec<(usize, usize)> = (0..ncols)
        .flat_map(|i| ((i + 1)..ncols).map(move |j| (i, j)))
        .collect();

    let results: Vec<(usize, usize, f64)> = pairs
        .par_iter()
        .map(|&(i, j)| {
            let mut mismatches = 0;
            for k in 0..nrows {
                if mat.get(k, i) != mat.get(k, j) {
                    mismatches += 1;
                }
            }
            let dist = mismatches as f64 / nrows as f64;
            (i, j, dist)
        })
        .collect();

    for (i, j, dist) in results {
        res[(i, j)] = dist;
        res[(j, i)] = dist;
    }

    res
}

/////////////////////////////////
// Topological overlap measure //
/////////////////////////////////

/// Enum for the TOM function
#[derive(Debug)]
pub enum TomType {
    Version1,
    Version2,
}

/// Parsing the TOM type
pub fn parse_tom_types(s: &str) -> Option<TomType> {
    match s.to_lowercase().as_str() {
        "v1" => Some(TomType::Version1),
        "v2" => Some(TomType::Version2),
        _ => None,
    }
}

/// Calculates the topological overlap measure (TOM) for a given affinity matrix
///
/// The TOM quantifies the relative interconnectedness of two nodes in a network by measuring
/// how much they share neighbors relative to their connectivity. Higher TOM values indicate
/// nodes that are part of the same module or cluster.
///
/// ### Params
///
/// * `affinity_mat` - Symmetric affinity/adjacency matrix
/// * `signed` - Whether to use signed (absolute values) or unsigned connectivity
/// * `tom_type` - Algorithm version (Version1 or Version2)
///
/// ### Returns
/// Symmetric TOM matrix with values in [0,1] representing topological overlap
///
/// ### Mathematical Formulation
///
/// #### Connectivity
/// For node i: k_i = Σ_j |a_ij| (signed) or k_i = Σ_j a_ij (unsigned)
///
/// #### Shared Neighbors
/// For nodes i,j: l_ij = Σ_k (a_ik * a_kj) - a_ii*a_ij - a_ij*a_jj
///
/// #### TOM Calculation
///
/// **Version 1:**
/// - Numerator: a_ij + l_ij
/// - Denominator (unsigned): min(k_i, k_j) + 1 - a_ij
/// - Denominator (signed): min(k_i, k_j) + 1 - a_ij (if a_ij ≥ 0) or min(k_i, k_j) + 1 + a_ij (if a_ij < 0)
/// - TOM_ij = numerator / denominator
///
/// **Version 2:**
/// - Divisor (unsigned): min(k_i, k_j) + a_ij
/// - Divisor (signed): min(k_i, k_j) + a_ij (if a_ij ≥ 0) or min(k_i, k_j) - a_ij (if a_ij < 0)
/// - TOM_ij = 0.5 * (a_ij + l_ij/divisor)
pub fn calc_tom(affinity_mat: MatRef<f64>, signed: bool, tom_type: TomType) -> Mat<f64> {
    let n = affinity_mat.nrows();
    let mut tom_mat = Mat::<f64>::zeros(n, n);

    let connectivity = if signed {
        (0..n)
            .map(|i| (0..n).map(|j| affinity_mat.get(i, j).abs()).sum())
            .collect::<Vec<f64>>()
    } else {
        col_sums(affinity_mat.as_ref())
    };

    // Pre-compute for speed-ups -> Massive difference and good call @Claude
    let dot_products = affinity_mat.as_ref() * affinity_mat.as_ref();

    for i in 0..n {
        // Only upper triangle for speed...
        for j in (i + 1)..n {
            let a_ij = affinity_mat.get(i, j);

            // shared_neighbors = dot_product - excluded terms (k=i and k=j)
            let shared_neighbours = dot_products.get(i, j)
                - affinity_mat.get(i, i) * affinity_mat.get(i, j)
                - affinity_mat.get(i, j) * affinity_mat.get(j, j);

            let f_ki_kj = connectivity[i].min(connectivity[j]);

            let tom_value = match tom_type {
                TomType::Version1 => {
                    let numerator = a_ij + shared_neighbours;
                    let denominator = if signed {
                        if *a_ij >= 0.0 {
                            f_ki_kj + 1.0 - a_ij
                        } else {
                            f_ki_kj + 1.0 + a_ij
                        }
                    } else {
                        f_ki_kj + 1.0 - a_ij
                    };
                    numerator / denominator
                }
                TomType::Version2 => {
                    let divisor = if signed {
                        if *a_ij >= 0.0 {
                            f_ki_kj + a_ij
                        } else {
                            f_ki_kj - a_ij
                        }
                    } else {
                        f_ki_kj + a_ij
                    };
                    let neighbours = shared_neighbours / divisor;
                    0.5 * (a_ij + neighbours)
                }
            };

            tom_mat[(i, j)] = tom_value;
            tom_mat[(j, i)] = tom_value;
        }
    }

    tom_mat
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
