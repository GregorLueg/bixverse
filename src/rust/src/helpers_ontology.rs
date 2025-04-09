use extendr_api::prelude::*;
use std::collections::{HashMap, HashSet};

/// Enum for the OntoSim types
#[derive(Debug)]
pub enum OntoSimType {
    Resnik,
    Lin,
}

/// Parsing the OntoSim types
pub fn parse_onto_sim(s: &str) -> Option<OntoSimType> {
    match s.to_lowercase().as_str() {
        "resnik" => Some(OntoSimType::Resnik),
        "lin" => Some(OntoSimType::Lin),
        _ => None,
    }
}

/// Get the most informative common ancestor
fn get_mica(
    t1: &str,
    t2: &str,
    ancestor_map: &HashMap<String, HashSet<String>>,
    info_content_map: &HashMap<String, f64>,
) -> f64 {
    // Default Hashmap to avoid all types of tests here...
    let default_hash: HashSet<String> =
        HashSet::from_iter(std::iter::once("I have no ancestors".to_string()));
    let ancestor_1 = ancestor_map.get(t1).unwrap_or(&default_hash);
    let ancestor_2 = ancestor_map.get(t2).unwrap_or(&default_hash);
    ancestor_1
        .intersection(ancestor_2)
        .map(|ancestor| info_content_map.get(ancestor).cloned().unwrap_or(0.0))
        .fold(0.0, f64::max)
}

/// Calculate the Lin similarity
pub fn lin_similarity(
    t1: &str,
    t2: &str,
    ancestor_map: &HashMap<String, HashSet<String>>,
    info_content_map: &HashMap<String, f64>,
) -> f64 {
    let mica = get_mica(t1, t2, ancestor_map, info_content_map);
    let t1_ic = info_content_map.get(t1).unwrap_or(&1.0);
    let t2_ic = info_content_map.get(t2).unwrap_or(&1.0);
    let max_value = f64::max(*t1_ic, *t2_ic);
    mica / max_value
}

/// Calculate the Resnik similarity (normalised)
pub fn resnik_similarity(
    t1: &str,
    t2: &str,
    max_ic: &f64,
    ancestor_map: &HashMap<String, HashSet<String>>,
    info_content_map: &HashMap<String, f64>,
) -> f64 {
    let mica = get_mica(t1, t2, ancestor_map, info_content_map);
    mica / max_ic
}

/// Calculate the Resnik and Lin similarity in one go
pub fn resnik_and_lin_sim(
    t1: &str,
    t2: &str,
    max_ic: &f64,
    ancestor_map: &HashMap<String, HashSet<String>>,
    info_content_map: &HashMap<String, f64>,
) -> (f64, f64) {
    let mica = get_mica(t1, t2, ancestor_map, info_content_map);
    let t1_ic = info_content_map.get(t1).unwrap_or(&1.0);
    let t2_ic = info_content_map.get(t2).unwrap_or(&1.0);
    let max_value = f64::max(*t1_ic, *t2_ic);
    let lin_sim = mica / max_value;
    let resnik_sim = mica / max_ic;
    (resnik_sim, lin_sim)
}

/// Transform an R list that hopefully contains the IC into a HashMap
/// of floats
pub fn ic_list_to_ic_hashmap(r_list: List) -> HashMap<String, f64> {
    let mut hashmap = HashMap::with_capacity(r_list.len());
    for (name, x) in r_list {
        let name = name.to_string();
        let ic_val = x.as_real().unwrap_or(0.0);
        hashmap.insert(name, ic_val);
    }
    hashmap
}
