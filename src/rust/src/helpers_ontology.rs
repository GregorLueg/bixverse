use crate::utils_rust::flatten_vector;
use extendr_api::prelude::*;
use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::Graph;
use rand::prelude::*;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

///////////////////////////
// Semantic similarities //
///////////////////////////

/// Structure to store the Ontology similarity results
#[derive(Clone, Debug)]
pub struct OntoSimRes<'a> {
    pub t1: &'a str,
    pub t2: &'a str,
    pub sim: f64,
}

/// Enum to store the different similarity types
#[derive(Clone, Debug)]
pub enum OntoSimType {
    Resnik,
    Lin,
    Combined,
}

/// Parsing the Onto Similarity types
fn parse_onto_sim_type(s: &str) -> Option<OntoSimType> {
    match s.to_lowercase().as_str() {
        "resnik" => Some(OntoSimType::Resnik),
        "lin" => Some(OntoSimType::Lin),
        "combined" => Some(OntoSimType::Combined),
        _ => None,
    }
}

/// Get the most informative common ancestors
fn get_mica(
    t1: &str,
    t2: &str,
    ancestor_map: &FxHashMap<String, FxHashSet<String>>,
    info_content_map: &FxHashMap<String, f64>,
) -> f64 {
    let default_hash: FxHashSet<String> =
        FxHashSet::from_iter(std::iter::once("I have no ancestors".to_string()));
    let ancestor_1 = ancestor_map.get(t1).unwrap_or(&default_hash);
    let ancestor_2 = ancestor_map.get(t2).unwrap_or(&default_hash);
    let mica = ancestor_1
        .intersection(ancestor_2)
        .map(|ancestor| info_content_map.get(ancestor).unwrap())
        .copied()
        .fold(0.0, f64::max);
    mica
}

/// Calculate the Resnik semantic similarity
fn calculate_resnik<'a>(
    t1: &'a str,
    t2: &'a str,
    ancestor_map: &FxHashMap<String, FxHashSet<String>>,
    info_content_map: &FxHashMap<String, f64>,
) -> OntoSimRes<'a> {
    let sim = get_mica(t1, t2, ancestor_map, info_content_map);
    OntoSimRes { t1, t2, sim }
}

/// Calculate the Lin semantic similarity
fn calculate_lin<'a>(
    t1: &'a str,
    t2: &'a str,
    ancestor_map: &FxHashMap<String, FxHashSet<String>>,
    info_content_map: &FxHashMap<String, f64>,
) -> OntoSimRes<'a> {
    let mica = get_mica(t1, t2, ancestor_map, info_content_map);
    let t1_ic = info_content_map.get(t1).unwrap_or(&1.0);
    let t2_ic = info_content_map.get(t2).unwrap_or(&1.0);
    let sim = 2.0 * mica / (t1_ic + t2_ic);
    OntoSimRes { t1, t2, sim }
}

/// Calculate for the combined Resnik/Lin similarity
fn calculate_combined_sim<'a>(
    t1: &'a str,
    t2: &'a str,
    max_ic: &f64,
    ancestor_map: &FxHashMap<String, FxHashSet<String>>,
    info_content_map: &FxHashMap<String, f64>,
) -> OntoSimRes<'a> {
    let mica = get_mica(t1, t2, ancestor_map, info_content_map);
    let t1_ic = info_content_map.get(t1).unwrap_or(&1.0);
    let t2_ic = info_content_map.get(t2).unwrap_or(&1.0);
    let lin_sim = 2.0 * mica / (t1_ic + t2_ic);
    let resnik_sim = mica / max_ic;
    let sim = (lin_sim + resnik_sim) / 2.0;
    OntoSimRes { t1, t2, sim }
}

/// Calculate the semantic similarity given two terms, the specified similarity
/// type and other needed information.
pub fn get_single_onto_sim<'a>(
    t1: &'a str,
    t2: &'a str,
    sim_type: &str,
    max_ic: &f64,
    ancestor_map: &FxHashMap<String, FxHashSet<String>>,
    info_content_map: &FxHashMap<String, f64>,
) -> Result<OntoSimRes<'a>> {
    let onto_sim_type = parse_onto_sim_type(sim_type)
        .ok_or_else(|| format!("Invalid Ontology Similarity Type: {}", sim_type))?;

    let res = match onto_sim_type {
        OntoSimType::Resnik => calculate_resnik(t1, t2, ancestor_map, info_content_map),
        OntoSimType::Lin => calculate_lin(t1, t2, ancestor_map, info_content_map),
        OntoSimType::Combined => {
            calculate_combined_sim(t1, t2, max_ic, ancestor_map, info_content_map)
        }
    };

    Ok(res)
}

/// Calculate the ontological similarity in an efficient manner for a set of terms
pub fn calculate_onto_sim<'a>(
    terms_split: &'a Vec<(String, &[String])>,
    sim_type: &str,
    ancestors_map: FxHashMap<String, FxHashSet<String>>,
    ic_map: FxHashMap<String, f64>,
) -> Vec<OntoSimRes<'a>> {
    let max_ic = ic_map
        .values()
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();

    let onto_sim: Vec<Vec<OntoSimRes<'_>>> = terms_split
        .par_iter()
        .map(|(t1, others)| {
            let mut sim_vec: Vec<OntoSimRes<'_>> = Vec::with_capacity(others.len());
            others.iter().for_each(|t2| {
                let sim_res =
                    get_single_onto_sim(t1, t2, sim_type, max_ic, &ancestors_map, &ic_map);
                sim_vec.push(sim_res.unwrap());
            });

            sim_vec
        })
        .collect();

    flatten_vector(onto_sim)
}

/// Transform an R list that hopefully contains the IC into a HashMap of floats
pub fn ic_list_to_ic_hashmap(r_list: List) -> FxHashMap<String, f64> {
    let mut hashmap = FxHashMap::default();
    for (name, x) in r_list {
        let name = name.to_string();
        let ic_val = x.as_real().unwrap_or(0.0);
        hashmap.insert(name, ic_val);
    }
    hashmap
}

/// Calculates the critical value
pub fn calculate_critval(values: &[f64], sample_size: usize, alpha: &f64, seed: usize) -> f64 {
    let mut rng = StdRng::seed_from_u64(seed as u64);
    let mut random_sample: Vec<f64> = (0..sample_size)
        .map(|_| {
            let index = rng.random_range(0..values.len());
            values[index]
        })
        .collect();
    random_sample.sort_by(|a, b| b.partial_cmp(a).unwrap());
    let index = (alpha * random_sample.len() as f64).ceil() as usize;
    random_sample[index + 1]
}

///////////////////////
// DAG-based methods //
///////////////////////
