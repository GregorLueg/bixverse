use extendr_api::prelude::*;
use faer::{unzip, zip, Mat};
use rustc_hash::FxHashSet;
use std::collections::BTreeMap;

use crate::helpers::linalg::cor;
use crate::utils::general::*;
use crate::utils::r_rust_interface::NamedMatrix;
use crate::utils::utils_stats::set_similarity;

////////////////////////
// Results structures //
////////////////////////

/// Structure for an Rbh triplet Result
///
/// ### Fields
///
/// * `t1` - Name of term 1 of the ontology
/// * `t2` - Name of term 2 of the ontology
/// * `sim` - Calculated similarity
#[derive(Clone, Debug)]
pub struct RbhTripletStruc<'a> {
    pub t1: &'a str,
    pub t2: &'a str,
    pub sim: f64,
}

/// Structure to store the RBH results.
///
/// ### Fields
///
/// * `origin` - Name of the origin data set
/// * `target` - Name of the target data set
/// * `origin_modules` - Names of the origin modules/gene sets
/// * `target_modules` - Names of the target modules/gene sets
/// * `similarities` - Similarities between the modules/gene sets
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct RbhResult {
    pub origin: String,
    pub target: String,
    pub origin_modules: Vec<String>,
    pub target_modules: Vec<String>,
    pub similarities: Vec<f64>,
}

////////////////////
// Set similarity //
////////////////////

/// Calculates the reciprocal best hits based on set similarities.
///
/// Function will calculate the set similarities (Jaccard or Overlap coefficient)
/// between all of the gene sets between the two data sets and calculate the
/// reciprocal best hits based on this similarity matrix.
///
/// ### Params
///
/// * `origin_modules` - A BTreeMap containing the identified modules of the
///   the origin data set.
/// * `target_modules` - A BTreeMap containing the identified modules of the
///   the target data set.
/// * `overlap_coefficient` - Shall the overlap coefficient be used instead of
///   Jaccard similarity.
/// * `min_similarity` - Minimum similarity to be returned
///
/// ### Returns
///
/// A vector of `RbhTripletStruc`.
pub fn calculate_rbh_set<'a>(
    origin_modules: &'a BTreeMap<String, FxHashSet<String>>,
    target_modules: &'a BTreeMap<String, FxHashSet<String>>,
    overlap_coefficient: bool,
    min_similarity: f64,
) -> Vec<RbhTripletStruc<'a>> {
    let names_targets: Vec<&String> = target_modules.keys().collect();
    let names_origin: Vec<&String> = origin_modules.keys().collect();

    let similarities_flat: Vec<Vec<f64>> = origin_modules
        .values()
        .map(|v1| {
            let v1_refs: FxHashSet<&String> = v1.iter().collect();
            target_modules
                .values()
                .map(|v2| {
                    let v2_refs: FxHashSet<&String> = v2.iter().collect();
                    set_similarity(&v1_refs, &v2_refs, overlap_coefficient)
                })
                .collect()
        })
        .collect();

    let mat_data: Vec<f64> = flatten_vector(similarities_flat);

    let max_sim = array_max(&mat_data);

    let result = if max_sim < min_similarity {
        vec![RbhTripletStruc {
            t1: "NA",
            t2: "NA",
            sim: 0.0,
        }]
    } else {
        let nrow = names_origin.len();
        let ncol = names_targets.len();

        let sim_mat = Mat::from_fn(nrow, ncol, |i, j| mat_data[j + i * ncol]);

        let row_maxima: Vec<&f64> = sim_mat
            .row_iter()
            .map(|x| {
                let row: Vec<&f64> = x.iter().collect();
                array_max(&row)
            })
            .collect();

        let col_maxima: Vec<&f64> = sim_mat
            .col_iter()
            .map(|x| {
                let col: Vec<&f64> = x.iter().collect();
                array_max(&col)
            })
            .collect();

        let mut matching_pairs: Vec<RbhTripletStruc> = Vec::new();

        for r in 0..nrow {
            for c in 0..ncol {
                let value = sim_mat[(r, c)];
                if value == *row_maxima[r] && value == *col_maxima[c] {
                    let triplet = RbhTripletStruc {
                        t1: names_origin[r],
                        t2: names_targets[c],
                        sim: value,
                    };

                    matching_pairs.push(triplet)
                }
            }
        }

        if !matching_pairs.is_empty() {
            matching_pairs
        } else {
            vec![RbhTripletStruc {
                t1: "NA",
                t2: "NA",
                sim: 0.0,
            }]
        }
    };

    result
}

///////////////////////
// Correlation based //
///////////////////////

/// Transforms a list of R matrices into a vector of R matrices
///
/// ### Params
///
/// * `matrix_list` - R List of matrices
///
/// ### Returns
///
/// A vector of tuples with the name of the list element and the R matrix.
pub fn r_matrix_list_to_vec(matrix_list: List) -> Vec<(String, RArray<f64, [usize; 2]>)> {
    matrix_list
        .iter()
        .map(|(n, obj)| (n.to_string(), obj.as_matrix().unwrap()))
        .collect()
}

/// Take a vector of R matrices and generate a BTreeMap of NamedMatrices
///
/// ### Params
///
/// * `matrix_vector` - Slice of tuples with the first element representing the name
///   and the second the R matrix
///
/// ### Returns
///
/// A BTreeMap of `NamedMatrix` objects.
pub fn r_matrix_vec_to_btree_list(
    matrix_vector: &[(String, RArray<f64, [usize; 2]>)],
) -> BTreeMap<String, NamedMatrix<'_>> {
    let mut result = BTreeMap::new();
    for (name, matrix) in matrix_vector {
        let named_mat = NamedMatrix::new(matrix);
        result.insert(name.clone(), named_mat);
    }

    result
}

/// Calculate the RBH based on correlation of two NamedMatrices
///
/// The function will intersect into shared features and calculate the correlation
/// matrix and subsequently reciprocal best hits based on the absolute correlation.
///
/// ### Params
///
/// * `x1` - `NamedMatrix` of the origin data
/// * `x2` - `NamedMatrix` of the target data
/// * `spearman` - Shall Spearman correlations be used.
///
/// ### Returns
///
/// A vector of `RbhTripletStruc`.
pub fn calculate_rbh_cor<'a>(
    x1: &'a NamedMatrix<'a>,
    x2: &'a NamedMatrix<'a>,
    spearman: bool,
) -> Vec<RbhTripletStruc<'a>> {
    let row_names_1: FxHashSet<String> = x1.row_names.keys().cloned().collect();
    let row_names_2: FxHashSet<String> = x2.row_names.keys().cloned().collect();

    // Now these references will live as long as 'a because they reference x1 and x2
    let names_targets: Vec<&str> = x1.col_names.keys().map(|s| s.as_str()).collect();
    let names_origin: Vec<&str> = x2.col_names.keys().map(|s| s.as_str()).collect();

    let intersecting_rows: Vec<String> = row_names_1.intersection(&row_names_2).cloned().collect();

    // Early return if there are no intersecting rows
    if intersecting_rows.is_empty() {
        return vec![RbhTripletStruc {
            t1: "NA",
            t2: "NA",
            sim: 0.0,
        }];
    }

    let row_refs: Vec<&str> = intersecting_rows.iter().map(|s| s.as_str()).collect();
    let sub_x1 = x1.get_rows(&row_refs).unwrap();
    let sub_x2 = x2.get_rows(&row_refs).unwrap();

    let mut correlations = cor(&sub_x1.as_ref(), &sub_x2.as_ref(), spearman);

    zip!(correlations.as_mut()).for_each(|unzip!(x)| *x = x.abs());

    let row_maxima: Vec<&f64> = correlations
        .row_iter()
        .map(|x| {
            let row: Vec<&f64> = x.iter().collect();
            array_max(&row)
        })
        .collect();

    let col_maxima: Vec<&f64> = correlations
        .col_iter()
        .map(|x| {
            let col: Vec<&f64> = x.iter().collect();
            array_max(&col)
        })
        .collect();

    let mut matching_pairs: Vec<RbhTripletStruc<'a>> = Vec::new();
    let nrow = names_targets.len();
    let ncol = names_origin.len();

    for r in 0..nrow {
        for c in 0..ncol {
            let value = correlations[(r, c)];
            if value == *row_maxima[r] && value == *col_maxima[c] {
                let triplet = RbhTripletStruc {
                    t1: names_targets[r],
                    t2: names_origin[c],
                    sim: value,
                };
                matching_pairs.push(triplet);
            }
        }
    }

    if !matching_pairs.is_empty() {
        matching_pairs
    } else {
        vec![RbhTripletStruc {
            t1: "NA",
            t2: "NA",
            sim: 0.0,
        }]
    }
}
