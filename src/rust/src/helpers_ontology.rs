use crate::utils_rust::flatten_vector;

use extendr_api::prelude::*;
use faer::Mat;
use petgraph::graph::{DiGraph, NodeIndex};
use rand::prelude::*;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::VecDeque;
use std::sync::Mutex;

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

#[derive(Clone, Debug)]
pub struct FastOntology {
    term_to_idx: FxHashMap<String, NodeIndex>, // term id to node index
    idx_to_term: Vec<String>,
    graph: DiGraph<String, ()>,           // directed graph: parent > child
    ancestors: Vec<FxHashSet<NodeIndex>>, // ancestors
    topo_order: Vec<NodeIndex>,           // pre-computing topological order
}

impl FastOntology {
    /// Create a new ontology object from parents and child
    pub fn new(parents: &[String], children: &[String]) -> Result<Self> {
        if parents.len() != children.len() {
            return Err("Parent and child vectors must have same length".into());
        }

        let mut graph = DiGraph::new();
        let mut term_to_idx = FxHashMap::default();
        let mut idx_to_term = Vec::new();

        let mut all_terms = FxHashSet::default();
        all_terms.extend(parents.iter().cloned());
        all_terms.extend(children.iter().cloned());

        for term in all_terms {
            let idx = graph.add_node(term.clone());
            term_to_idx.insert(term.clone(), idx);
            idx_to_term.push(term);
        }

        for (parent, child) in parents.iter().zip(children.iter()) {
            let parent_idx = *term_to_idx.get(parent).unwrap();
            let child_idx = *term_to_idx.get(child).unwrap();
            graph.add_edge(parent_idx, child_idx, ());
        }

        let ancestors = Self::compute_ancestors(&graph, &term_to_idx);

        let topo_order = Self::compute_topological_order(&graph);

        Ok(FastOntology {
            term_to_idx,
            idx_to_term,
            graph,
            ancestors,
            topo_order,
        })
    }

    /// Calculate Wang similarity between two terms
    pub fn wang_sim(&self, term1: &str, term2: &str, w: f64) -> Option<f64> {
        let term_idx_1 = *self.term_to_idx.get(term1)?;
        let term_idx_2 = *self.term_to_idx.get(term2)?;

        if term_idx_1 == term_idx_2 {
            return Some(1.0);
        }

        let s_val_1 = self.calculate_s_values(term_idx_1, w);
        let s_val_2 = self.calculate_s_values(term_idx_2, w);

        let dag1_nodes = &self.ancestors[term_idx_1.index()];
        let dag2_nodes = &self.ancestors[term_idx_2.index()];

        let common_nodes: Vec<NodeIndex> = dag1_nodes.intersection(dag2_nodes).cloned().collect();

        if common_nodes.is_empty() {
            return Some(0.0);
        }

        let sv1: f64 = s_val_1.values().sum();
        let sv2: f64 = s_val_2.values().sum();

        let numerator: f64 = common_nodes
            .iter()
            .map(|&node_idx| {
                s_val_1.get(&node_idx).unwrap_or(&0.0) + s_val_2.get(&node_idx).unwrap_or(&0.0)
            })
            .sum();

        let denominator = sv1 + sv2;

        if denominator > 0.0 {
            Some(numerator / denominator)
        } else {
            Some(0.0)
        }
    }

    /// Calculate the similarity matrix
    pub fn calc_sim_matrix(&self, w: f64) -> Mat<f64> {
        let n = self.idx_to_term.len();
        let mut matrix: Mat<f64> = Mat::zeros(n, n);

        let matrix_mutex = Mutex::new(&mut matrix);

        (0..n).into_par_iter().for_each(|i| {
            let mut row_sim = Vec::with_capacity(n - i);
            for j in i..n {
                let sim = if i == j {
                    1.0
                } else {
                    let term1 = &self.idx_to_term[i];
                    let term2 = &self.idx_to_term[j];
                    self.wang_sim(term1, term2, w).unwrap_or(0.0)
                };
                row_sim.push((j, sim));
            }

            {
                let mut matrix_guard = matrix_mutex.lock().unwrap();
                for (j, sim) in row_sim {
                    matrix_guard[(i, j)] = sim;
                    if i != j {
                        matrix_guard[(j, i)] = sim;
                    }
                }
            }
        });

        matrix
    }

    /// Pre-compute the ancestors via a BFS
    fn compute_ancestors(
        graph: &DiGraph<String, ()>,
        term_to_idx: &FxHashMap<String, NodeIndex>,
    ) -> Vec<FxHashSet<NodeIndex>> {
        let mut ancestors = vec![FxHashSet::default(); graph.node_count()];

        for (_, &node_idx) in term_to_idx.iter() {
            let mut visited = FxHashSet::default();
            let mut queue = VecDeque::new();

            queue.push_back(node_idx);
            visited.insert(node_idx);

            while let Some(current) = queue.pop_front() {
                for parent_idx in graph.neighbors_directed(current, petgraph::Incoming) {
                    if !visited.contains(&parent_idx) {
                        visited.insert(parent_idx);
                        queue.push_back(parent_idx);
                    }
                }
            }

            ancestors[node_idx.index()] = visited;
        }

        ancestors
    }

    /// Compute topological order
    fn compute_topological_order(graph: &DiGraph<String, ()>) -> Vec<NodeIndex> {
        petgraph::algo::toposort(graph, None)
            .unwrap_or_else(|_| panic!("Ontology contains cycles"))
            .into_iter()
            .rev()
            .collect()
    }

    /// Calculate the S-values for a specific term's DAG
    fn calculate_s_values(&self, term_idx: NodeIndex, w: f64) -> FxHashMap<NodeIndex, f64> {
        let dag_nodes = &self.ancestors[term_idx.index()];
        let mut s_values = FxHashMap::default();

        s_values.insert(term_idx, 1.0);

        for &node_idx in &self.topo_order {
            if !dag_nodes.contains(&node_idx) || node_idx == term_idx {
                continue;
            }

            let mut max_contributions: f64 = 0.0;
            for child_idx in self.graph.neighbors_directed(node_idx, petgraph::Outgoing) {
                if dag_nodes.contains(&child_idx) {
                    if let Some(&child_s_value) = s_values.get(&child_idx) {
                        max_contributions = max_contributions.max(w * child_s_value);
                    }
                }
            }

            if max_contributions > 0.0 {
                s_values.insert(node_idx, max_contributions);
            }
        }

        s_values
    }
}
