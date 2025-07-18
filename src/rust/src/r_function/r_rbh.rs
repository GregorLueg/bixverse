use extendr_api::prelude::*;

use rayon::prelude::*;

use crate::helpers::rbh::*;
use crate::utils::general::flatten_vector;
use crate::utils::r_rust_interface::{r_nested_list_to_btree_nest, NestedBtreeMap};

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
pub struct RbhResult {
    pub origin: String,
    pub target: String,
    pub origin_modules: Vec<String>,
    pub target_modules: Vec<String>,
    pub similarities: Vec<f64>,
}

/// Generate reciprocal best hits based on set similarities
///
/// @description This function takes a nested list that contains gene modules/
/// sets derived from various methods and generate identifies reciprocal best
/// hits between gene modules/sets across the different origins. WARNING!
/// Incorrect use can cause kernel crashes. Wrapper around the Rust functions
/// with type checks are provided in the package.
///
/// @param module_list A nested named list. The outer list should contain the
/// origin of the gene modules, the inner list the names of the gene modules and
/// the respective genes in them.
/// @param overlap_coefficient Shall the overlap coefficient instead of the
/// Jaccard similarity be used.
/// @param min_similarity Minimum similarity that should exist between any two
/// given gene modules to actually calculate RBH pairs.
///
/// @return A list containing:
///  \itemize{
///   \item origin - The name of the origin of the gene modules.
///   \item target - The name of the target of the gene modules.
///   \item comparisons - Integer vector indicating how many RBH hits were
///   identified in this comparison
///   \item origin_modules - Names of the gene modules from the origin.
///   \item target_modules - Names of the gene modules from the target.
///   \item similarity - The similarities between the two respective gene
///   modules.
/// }
///
/// @export
#[extendr]
fn rs_rbh_sets(
    module_list: List,
    overlap_coefficient: bool,
    min_similarity: f64,
) -> extendr_api::Result<List> {
    let module_list: NestedBtreeMap = r_nested_list_to_btree_nest(module_list)?;
    // Pull out all the keys
    let origins: Vec<String> = module_list.keys().cloned().collect();

    let origins_split: Vec<(String, &[String])> = origins
        .iter()
        .enumerate()
        .map(|(i, first)| (first.clone(), &origins[i + 1..]))
        .take_while(|(_, rest)| !rest.is_empty())
        .collect();

    let rbh_results: Vec<Vec<RbhResult>> = origins_split
        .par_iter()
        .map(|(origin_module, target_modules)| {
            // Parallel iterator starts here
            let origin_module_data = module_list.get(origin_module).unwrap();

            // Iterate over the remaining target modules
            let res: Vec<RbhResult> = target_modules
                .iter()
                .map(|target| {
                    let target_module_data = module_list.get(target).unwrap();

                    let rbh_res = calculate_rbh_set(
                        origin_module_data,
                        target_module_data,
                        overlap_coefficient,
                        min_similarity,
                    );

                    let mut origin_modules = Vec::new();
                    let mut target_modules = Vec::new();
                    let mut similarities = Vec::new();

                    for res in rbh_res {
                        if res.sim > min_similarity {
                            origin_modules.push(res.t1.to_string());
                            target_modules.push(res.t2.to_string());
                            similarities.push(res.sim);
                        }
                    }

                    RbhResult {
                        origin: origin_module.to_string(),
                        target: target.to_string(),
                        origin_modules,
                        target_modules,
                        similarities,
                    }
                })
                .collect();

            res
        })
        .collect();

    // Flatten and extract relevant data.
    let rbh_results_flatten: Vec<_> = flatten_vector(rbh_results);

    let mut origin = Vec::new();
    let mut target = Vec::new();
    let mut comparisons = Vec::new();
    let mut origin_modules = Vec::new();
    let mut target_modules = Vec::new();
    let mut similarity = Vec::new();

    for module in rbh_results_flatten {
        origin.push(module.origin);
        target.push(module.target);
        origin_modules.push(module.origin_modules);
        target_modules.push(module.target_modules);
        similarity.push(module.similarities.clone());
        comparisons.push(module.similarities.len());
    }

    let origin_modules: Vec<_> = flatten_vector(origin_modules);
    let target_modules: Vec<_> = flatten_vector(target_modules);
    let similarity: Vec<_> = flatten_vector(similarity);

    Ok(list!(
        origin = origin,
        target = target,
        comparisons = comparisons,
        origin_modules = origin_modules,
        target_modules = target_modules,
        similarity = similarity
    ))
}

/// Generate reciprocal best hits based on correlations
///
/// @description This function takes list of (named) matrices which represent
/// for example matrix factorisation results you wish to identify reciprocal
/// best hits (RBH) for. The rows need to represent the features and the columns
/// the parts you wish to calculate the RBH for.
///
/// @param module_matrices A list of named matrices. Rows represent features
/// and columns the samples you wish to calculate the correlations for.
/// @param spearman Shall Spearman correlation be used.
/// @param min_similarity Minimum (absolute) correlations that needs to exist
/// between two terms.
///
/// @return A list containing:
///  \itemize{
///   \item origin - The name of the origin of the gene modules.
///   \item target - The name of the target of the gene modules.
///   \item comparisons - Integer vector indicating how many RBH hits were
///   identified in this comparison
///   \item origin_modules - Names of the gene modules from the origin.
///   \item target_modules - Names of the gene modules from the target.
///   \item similarity - The similarities between the two respective gene
///   modules.
/// }
///
/// @export
#[extendr]
fn rs_rbh_cor(
    module_matrices: List,
    spearman: bool,
    min_similarity: f64,
) -> extendr_api::Result<List> {
    let matrix_vec = r_matrix_list_to_vec(module_matrices);

    let matrix_map = r_matrix_vec_to_btree_list(&matrix_vec);

    let origins: Vec<String> = matrix_map.keys().cloned().collect();

    let origins_split: Vec<(String, &[String])> = origins
        .iter()
        .enumerate()
        .map(|(i, first)| (first.clone(), &origins[i + 1..]))
        .take_while(|(_, rest)| !rest.is_empty())
        .collect();

    let rbh_results: Vec<Vec<RbhResult>> = origins_split
        .par_iter()
        .map(|(origin_module, target_modules)| {
            let origin_mat = matrix_map.get(origin_module).unwrap();
            let res: Vec<RbhResult> = target_modules
                .iter()
                .map(|target| {
                    let target_matrix = matrix_map.get(target).unwrap();
                    let result = calculate_rbh_cor(origin_mat, target_matrix, spearman);
                    let mut origin_modules: Vec<String> = Vec::new();
                    let mut target_modules: Vec<String> = Vec::new();
                    let mut similarities: Vec<f64> = Vec::new();

                    for res in result {
                        if res.sim > min_similarity {
                            origin_modules.push(res.t1.to_string());
                            target_modules.push(res.t2.to_string());
                            similarities.push(res.sim);
                        }
                    }

                    RbhResult {
                        origin: origin_module.to_string(),
                        target: target.to_string(),
                        origin_modules,
                        target_modules,
                        similarities,
                    }
                })
                .collect();

            res
        })
        .collect();

    // Flatten and extract relevant data.
    let rbh_results_flatten: Vec<_> = flatten_vector(rbh_results);

    let mut origin = Vec::new();
    let mut target = Vec::new();
    let mut comparisons = Vec::new();
    let mut origin_modules = Vec::new();
    let mut target_modules = Vec::new();
    let mut similarity = Vec::new();

    for module in rbh_results_flatten {
        origin.push(module.origin);
        target.push(module.target);
        origin_modules.push(module.origin_modules);
        target_modules.push(module.target_modules);
        similarity.push(module.similarities.clone());
        comparisons.push(module.similarities.len());
    }

    let origin_modules: Vec<_> = flatten_vector(origin_modules);
    let target_modules: Vec<_> = flatten_vector(target_modules);
    let similarity: Vec<_> = flatten_vector(similarity);

    Ok(list!(
        origin = origin,
        target = target,
        comparisons = comparisons,
        origin_modules = origin_modules,
        target_modules = target_modules,
        similarity = similarity
    ))
}

extendr_module! {
    mod r_rbh;
    fn rs_rbh_sets;
    fn rs_rbh_cor;
}
