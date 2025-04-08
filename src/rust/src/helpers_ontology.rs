use extendr_api::prelude::*;
use std::collections::{HashMap, HashSet};

use crate::utils_rust::array_f64_max;

/// Enum for the ICA types
#[derive(Debug)]
pub enum OntoSimType {
    Resnik,
    Lin,
}

/// Parsing the ICA types
pub fn parse_onto_sim(s: &str) -> Option<OntoSimType> {
    match s.to_lowercase().as_str() {
        "resnik" => Some(OntoSimType::Resnik),
        "lin" => Some(OntoSimType::Lin),
        _ => None,
    }
}

/// Get the most informative common ancestor
pub fn get_mica(
    t1: &str,
    t2: &str,
    ancestor_map: &HashMap<String, HashSet<String>>,
    info_content_map: &HashMap<String, f64>,
) -> f64 {
    // Default Hashmap to avoid all types of tests here...
    let default_s = vec!["I have no ancestors".to_string()];
    let default_hash: HashSet<_> = default_s.into_iter().collect();
    let ancestor_1 = ancestor_map.get(t1).unwrap_or(&default_hash);
    let ancestor_2 = ancestor_map.get(t2).unwrap_or(&default_hash);
    let common_ancestors: Vec<String> = ancestor_1.intersection(ancestor_2).cloned().collect();
    if !common_ancestors.is_empty() {
        let mut info_content = Vec::new();
        for ancestor in common_ancestors {
            let ic = info_content_map.get(&ancestor).cloned().unwrap_or(0.0);
            info_content.push(ic);
        }
        array_f64_max(&info_content)
    } else {
        0.0
    }
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
    let mut hashmap = HashMap::new();
    for (name, x) in r_list {
        let name = name.to_string();
        let ic_val = x.as_real().unwrap_or(0.0);
        hashmap.insert(name, ic_val);
    }
    hashmap
}
