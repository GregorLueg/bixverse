use extendr_api::prelude::*;

use crate::helpers_ontology::*;
use crate::utils_r_rust::{faer_to_r_matrix, r_list_to_hashmap_set};
use crate::utils_rust::faer_mat_to_upper_triangle;

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
fn rs_onto_semantic_sim(
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

/// Calculate the Wang similarity for an ontology
///
/// @description This function calculates the Wang similarity for a given
/// ontology.
///
/// @param parents String vector. The names of the parents.
/// @param children String vector. The names of the childs. The length of
/// `parents` needs to be equal to `children`.
/// @param w Float. The w parameter for the ontology. Needs to be between
/// `0 < w < 1`.
/// @param flat_matrix Boolean. Shall only the upper triangle be returned.
///
/// @return A list with:
/// \itemize{
///   \item sim_mat - the Wang similarity matrix.
///   \item names - the row and column names for the calculated matrix.
/// }
///
/// @export
#[extendr]
fn rs_onto_sim_wang(
    parents: Vec<String>,
    children: Vec<String>,
    w: f64,
    flat_matrix: bool,
) -> extendr_api::Result<List> {
    let fast_onto = WangSimOntology::new(&parents, &children)?;

    let (sim_mat, names) = fast_onto.calc_sim_matrix(w);

    fast_onto.clear_cache();

    // Early return with flat matrix
    if flat_matrix {
        return Ok(list!(
            sim_mat = faer_mat_to_upper_triangle(sim_mat.as_ref(), 1),
            names = names
        ));
    };

    Ok(list!(
        sim_mat = faer_to_r_matrix(sim_mat.as_ref()),
        names = names
    ))
}

extendr_module! {
  mod fun_ontology;
  fn rs_onto_semantic_sim;
  fn rs_onto_sim_wang;
}
