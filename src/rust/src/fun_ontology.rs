use extendr_api::prelude::*;

use crate::helpers_ontology::*;
use crate::utils_r_rust::{faer_to_r_matrix, r_list_to_hashmap_set};

/// Calculate the semantic similarity in an ontology
///
/// @description This function calculates the specified semantic similarity and
/// returns the full vector (only calculating the upper triangle) for the given
/// similarity.
///
/// @param terms Vector of strings. The terms in the ontology you wish to screen.
/// @param sim_type String. Must be one of `c("resnik", "lin", "combined")`.
/// @param ancestor_list R list with names being the term and the elements in the
/// list the names of the ancestors.
/// @param ic_list R list with the names being the term and the elements the
/// information content of this given term. Needs to be a single float!
///
/// @return A vector containing all the desired similarity scores. This is
/// equivalent of the upper triangle of the similarity matrix.
///
/// @export
#[extendr]
fn rs_onto_similarity(
    terms: Vec<String>,
    sim_type: String,
    ancestor_list: List,
    ic_list: List,
) -> extendr_api::Result<Vec<f64>> {
    let ancestors_map = r_list_to_hashmap_set(ancestor_list)?;
    let ic_map = ic_list_to_ic_hashmap(ic_list);

    let terms_split: Vec<(String, &[String])> = terms
        .iter()
        .enumerate()
        .map(|(i, first)| (first.clone(), &terms[i + 1..]))
        .take_while(|(_, rest)| !rest.is_empty())
        .collect();

    let onto_sim = calculate_onto_sim(&terms_split, &sim_type, ancestors_map, ic_map);

    let mut final_sim = Vec::with_capacity(onto_sim.len());

    for sim_res in onto_sim.iter() {
        final_sim.push(sim_res.sim)
    }

    Ok(final_sim)
}

/// Calculate the semantic similarity in an ontology
///
/// @description This function calculates the specified semantic similarity and
/// returns the full vector (only calculating the upper triangle) for the given
/// similarity.
///
/// @param terms Vector of strings. The terms in the ontology you wish to screen.
/// @param sim_type String. Must be one of `c("resnik", "lin", "combined")`.
/// @param alpha Float. Must be between 0 to 1. The alpha parameter for calculating
/// the critival value.
/// @param ancestor_list R list with names being the term and the elements in the
/// list the names of the ancestors.
/// @param ic_list R list with the names being the term and the elements the
/// information content of this given term. Needs to be a single float!
/// @param iters Integer. Number of random iterations to use to estimate the
/// critical value.
/// @param seed Integer. Random seed for reproducibility purposes.
///
/// @return A list with:
/// \itemize{
///   \item term1 - Term 1
///   \item v - v matrix of the SVD.
///   \item s - Eigenvalues of the SVD.
/// }
///
/// @export
#[extendr]
fn rs_onto_similarity_filtered(
    terms: Vec<String>,
    sim_type: String,
    alpha: f64,
    ancestor_list: List,
    ic_list: List,
    iters: usize,
    seed: usize,
) -> extendr_api::Result<List> {
    let ancestors_map = r_list_to_hashmap_set(ancestor_list)?;
    let ic_map = ic_list_to_ic_hashmap(ic_list);

    let terms_split: Vec<(String, &[String])> = terms
        .iter()
        .enumerate()
        .map(|(i, first)| (first.clone(), &terms[i + 1..]))
        .take_while(|(_, rest)| !rest.is_empty())
        .collect();

    let onto_sim = calculate_onto_sim(&terms_split, &sim_type, ancestors_map, ic_map);

    let mut intermediate_sim = Vec::with_capacity(onto_sim.len());

    for sim_res in onto_sim.iter() {
        intermediate_sim.push(sim_res.sim)
    }

    let critval = calculate_critval(&intermediate_sim, iters, &alpha, seed);

    let mut term1 = Vec::new();
    let mut term2 = Vec::new();
    let mut final_sim = Vec::new();

    for sim_res in onto_sim.iter() {
        if sim_res.sim >= critval {
            term1.push(sim_res.t1.to_string());
            term2.push(sim_res.t2.to_string());
            final_sim.push(sim_res.sim);
        }
    }

    Ok(list!(
        term1 = term1,
        term2 = term2,
        filtered_sim = final_sim,
        critval = critval,
    ))
}

#[extendr]
fn rs_onto_sim_wang(
    parents: Vec<String>,
    children: Vec<String>,
    w: f64,
) -> extendr_api::Result<RArray<f64, [usize; 2]>> {
    let fast_onto = FastOntology::new(&parents, &children)?;

    let sim_mat = fast_onto.calc_sim_matrix(w);

    Ok(faer_to_r_matrix(sim_mat.as_ref()))
}

extendr_module! {
  mod fun_ontology;
  fn rs_onto_similarity;
  fn rs_onto_similarity_filtered;
  fn rs_onto_sim_wang;
}
