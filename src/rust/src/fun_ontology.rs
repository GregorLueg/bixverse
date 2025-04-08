use extendr_api::prelude::*;
use rayon::prelude::*;

use crate::helpers_ontology::*;
use crate::utils_r_rust::r_list_to_hashmap_set;
use crate::utils_rust::flatten_vector;

/// Calculate the semantic similarity in an ontology
///
/// @description This function calculates the Resnik or Lin similarity for a given ontology.
///
/// @param terms Vector of strings. The terms in the ontology you wish to screen.
/// @param ancestor_list R list with names being the term and the elements in the list the names
/// of the ancestors.
/// @param ic_list R list with the names being the term and the elements the information content
/// of this given term. Needs to be a single float!
/// @param similarity_type String. Need to be one of `c("resnik", "lin")`.
/// @param max_ic Double. The maximum information content observed in the data. This will return
/// normalised Resnik (i.e., scaled between 0 and 1). If set to 1, it will return the Resnik distance.
///
/// @return A list with:
/// \itemize{
///   \item terms - The supplied iterated terms.
///   \item similarities - Numeric vector of the similarities.
/// }
///
/// @export
#[extendr]
fn rs_onto_similarity(
    terms: Vec<String>,
    ancestor_list: List,
    ic_list: List,
    similarity_type: &str,
    max_ic: f64,
) -> extendr_api::Result<List> {
    let ancestors_map = r_list_to_hashmap_set(ancestor_list)?;
    let ic_map = ic_list_to_ic_hashmap(ic_list);

    let sim_type = parse_onto_sim(similarity_type)
        .ok_or_else(|| format!("Invalid Ontological similarity chosen: {}", similarity_type))?;

    let terms_split: Vec<(String, Vec<String>)> = terms
        .iter()
        .enumerate()
        .map(|(i, first)| {
            let rest: Vec<String> = terms.iter().skip(i + 1).map(|s| s.to_string()).collect();
            (first.to_string(), rest)
        })
        .take_while(|(_, rest)| !rest.is_empty())
        .collect();

    let onto_sim: Vec<Vec<f64>> = terms_split
        .par_iter()
        .map(|(t1, others)| {
            let sim_vec: Vec<f64> = others
                .par_iter()
                .map(|t2| match sim_type {
                    OntoSimType::Lin => lin_similarity(t1, t2, &ancestors_map, &ic_map),
                    OntoSimType::Resnik => {
                        resnik_similarity(t1, t2, &max_ic, &ancestors_map, &ic_map)
                    }
                })
                .collect();
            sim_vec
        })
        .collect();

    let onto_sim = flatten_vector(onto_sim);

    Ok(list!(terms = terms, similarities = onto_sim))
}

/// Calculate the Resnik and Lin semantic similarity
///
/// @description This function calculates the Resnik and Lin similarity for a given ontology in
/// one call.
///
/// @param terms Vector of strings. The terms in the ontology you wish to screen.
/// @param ancestor_list R list with names being the term and the elements in the list the names
/// of the ancestors.
/// @param ic_list R list with the names being the term and the elements the information content
/// of this given term. Needs to be a single float!
/// @param max_ic Double. The maximum information content observed in the data. This will return
/// normalised Resnik (i.e., scaled between 0 and 1). If set to 1, it will return the Resnik distance.
///
/// @return A list with:
/// \itemize{
///   \item terms - The supplied iterated terms.
///   \item resnik_sim - The Resnic similarities.
///   \item lin_sim - The Lin similarities.
/// }
///
/// @export
#[extendr]
fn rs_onto_similarity_both(
    terms: Vec<String>,
    ancestor_list: List,
    ic_list: List,
    max_ic: f64,
) -> extendr_api::Result<List> {
    let ancestors_map = r_list_to_hashmap_set(ancestor_list)?;
    let ic_map = ic_list_to_ic_hashmap(ic_list);

    let terms_split: Vec<(String, Vec<String>)> = terms
        .iter()
        .enumerate()
        .map(|(i, first)| {
            let rest: Vec<String> = terms.iter().skip(i + 1).map(|s| s.to_string()).collect();
            (first.to_string(), rest)
        })
        .take_while(|(_, rest)| !rest.is_empty())
        .collect();

    let onto_sim: Vec<Vec<(f64, f64)>> = terms_split
        .par_iter()
        .map(|(t1, others)| {
            let sim_vec: Vec<(f64, f64)> = others
                .par_iter()
                .map(|t2| resnik_and_lin_sim(t1, t2, &max_ic, &ancestors_map, &ic_map))
                .collect();
            sim_vec
        })
        .collect();

    let onto_sim = flatten_vector(onto_sim);

    let mut resnik_sim = Vec::new();
    let mut lin_sim = Vec::new();

    for (res, lin) in onto_sim.into_iter() {
        resnik_sim.push(res);
        lin_sim.push(lin);
    }

    Ok(list!(
        terms = terms,
        resnik_sim = resnik_sim,
        lin_sim = lin_sim
    ))
}

extendr_module! {
  mod fun_ontology;
  fn rs_onto_similarity;
  fn rs_onto_similarity_both;
}
