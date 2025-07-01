use crate::utils::general::flatten_vector;

use extendr_api::prelude::*;
use faer::Mat;
use petgraph::graph::{DiGraph, NodeIndex};
use rayon::prelude::*;
use rustc_hash::{FxBuildHasher, FxHashMap, FxHashSet};
use std::sync::{Arc, Mutex, RwLock};

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
    let mut hashmap = FxHashMap::with_capacity_and_hasher(r_list.len(), FxBuildHasher);
    for (name, x) in r_list {
        let name = name.to_string();
        let ic_val = x.as_real().unwrap_or(0.0);
        hashmap.insert(name, ic_val);
    }
    hashmap
}

///////////////////////
// DAG-based methods //
///////////////////////

/// Type alias for SValue Cache
/// This one needs the RwLock to be able to do the parallelisation on top
pub type SValueCache = RwLock<FxHashMap<(NodeIndex, u64), FxHashMap<NodeIndex, f64>>>;

/// Structure for calculating the Wang similarity on a given Ontology
#[derive(Debug)]
#[allow(dead_code)]
pub struct WangSimOntology {
    term_to_idx: FxHashMap<String, NodeIndex>, // term id to node index
    idx_to_term: Vec<String>,                  // the idx to term order (important)
    graph: DiGraph<String, ()>,                // directed graph: parent > child
    ancestors: Vec<FxHashSet<NodeIndex>>,      // ancestors (including self)
    topo_order: Vec<NodeIndex>,                // pre-computing topological order
    s_values_cache: SValueCache,               // pre-computed S-values
}

impl WangSimOntology {
    /// Create a new ontology object from parents and child
    /// Assumes two strings as input, one being the parents, the other being
    /// the children
    pub fn new(parents: &[String], children: &[String]) -> Result<Self> {
        if parents.len() != children.len() {
            return Err("Parent and child vectors must have same length".into());
        }

        let mut graph = DiGraph::new();
        let mut term_to_idx = FxHashMap::default();
        let mut idx_to_term = Vec::new();

        // More efficient collection of unique terms
        let mut all_terms = FxHashSet::with_capacity_and_hasher(
            (parents.len() + children.len()) * 2,
            FxBuildHasher,
        );
        all_terms.extend(parents.iter().cloned());
        all_terms.extend(children.iter().cloned());

        // Pre-allocate vectors with known capacity
        idx_to_term.reserve(all_terms.len());
        term_to_idx.reserve(all_terms.len());

        // Add nodes to the graph
        for term in all_terms {
            let idx = graph.add_node(term.clone());
            term_to_idx.insert(term.clone(), idx);
            idx_to_term.push(term);
        }

        // Add edges to the graph
        for (parent, child) in parents.iter().zip(children.iter()) {
            let parent_idx = *term_to_idx.get(parent).unwrap();
            let child_idx = *term_to_idx.get(child).unwrap();
            graph.add_edge(parent_idx, child_idx, ());
        }

        let topo_order = Self::compute_topological_order(&graph);
        let ancestors = Self::get_ancestors(&graph, &topo_order);

        Ok(WangSimOntology {
            term_to_idx,
            idx_to_term,
            graph,
            ancestors,
            topo_order,
            s_values_cache: RwLock::new(FxHashMap::default()),
        })
    }

    /// Calculate the similarity matrix with optimizations
    pub fn calc_sim_matrix(&self, w: f64) -> (Mat<f64>, Vec<String>) {
        let n = self.idx_to_term.len();
        let mut matrix: Mat<f64> = Mat::zeros(n, n);

        // Pre-compute all S-values once and store in Arc for safe sharing
        let w_key = w.to_bits();
        let all_s_values: Arc<Vec<FxHashMap<NodeIndex, f64>>> = Arc::new(
            (0..n)
                .into_par_iter()
                .map(|i| {
                    let term_idx = NodeIndex::new(i);
                    self.get_or_compute_s_values(term_idx, w, w_key)
                })
                .collect(),
        );

        // Process each row in parallel; mutex is a new one...
        let matrix_mutex = Mutex::new(&mut matrix);

        (0..n).into_par_iter().for_each(|i| {
            let mut row_values = Vec::with_capacity(n - i);

            // Calculate similarities for this row (only upper triangle)
            for j in i..n {
                let sim = if i == j {
                    1.0
                } else {
                    let term_idx_1 = NodeIndex::new(i);
                    let term_idx_2 = NodeIndex::new(j);
                    self.calculate_similarity_from_s_values(
                        term_idx_1,
                        term_idx_2,
                        &all_s_values[i],
                        &all_s_values[j],
                    )
                };
                row_values.push((j, sim));
            }

            // Write the computed values to the matrix
            {
                let mut matrix_guard = matrix_mutex.lock().unwrap();
                for (j, sim) in row_values {
                    matrix_guard[(i, j)] = sim;
                    if i != j {
                        matrix_guard[(j, i)] = sim; // Symmetric matrix
                    }
                }
            }
        });

        (matrix, self.idx_to_term.clone())
    }

    /// Clear the S-values cache (useful for memory management)
    pub fn clear_cache(&self) {
        let mut cache = self.s_values_cache.write().unwrap();
        cache.clear();
    }

    /// Get or compute S-values with caching
    fn get_or_compute_s_values(
        &self,
        term_idx: NodeIndex,
        w: f64,
        w_key: u64,
    ) -> FxHashMap<NodeIndex, f64> {
        let cache_key = (term_idx, w_key);

        // Try to read from cache first
        {
            let cache = self.s_values_cache.read().unwrap();
            if let Some(cached) = cache.get(&cache_key) {
                return cached.clone();
            }
        }

        // Compute if not in cache
        let s_values = self.calculate_s_values(term_idx, w);

        // Store in cache
        {
            let mut cache = self.s_values_cache.write().unwrap();
            cache.insert(cache_key, s_values.clone());
        }

        s_values
    }

    /// Calculate similarity from pre-computed S-values
    fn calculate_similarity_from_s_values(
        &self,
        term_idx_1: NodeIndex,
        term_idx_2: NodeIndex,
        s_val_1: &FxHashMap<NodeIndex, f64>,
        s_val_2: &FxHashMap<NodeIndex, f64>,
    ) -> f64 {
        let dag1_nodes = &self.ancestors[term_idx_1.index()];
        let dag2_nodes = &self.ancestors[term_idx_2.index()];

        // Start with smaller to make it faster
        let (smaller, larger) = if dag1_nodes.len() < dag2_nodes.len() {
            (dag1_nodes, dag2_nodes)
        } else {
            (dag2_nodes, dag1_nodes)
        };

        let common_nodes: Vec<NodeIndex> = smaller
            .iter()
            .filter(|node| larger.contains(node))
            .cloned()
            .collect();

        if common_nodes.is_empty() {
            return 0.0;
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
            numerator / denominator
        } else {
            0.0
        }
    }

    /// Get the ancestor terms of everything in the ontology
    fn get_ancestors(
        graph: &DiGraph<String, ()>,
        topo_order: &[NodeIndex],
    ) -> Vec<FxHashSet<NodeIndex>> {
        let mut ancestors = vec![FxHashSet::default(); graph.node_count()];

        // Process nodes in reverse topological order
        for &node_idx in topo_order.iter().rev() {
            let mut node_ancestors = FxHashSet::default();

            // Add self
            node_ancestors.insert(node_idx);
            for parent_idx in graph.neighbors_directed(node_idx, petgraph::Incoming) {
                node_ancestors.extend(&ancestors[parent_idx.index()]);
            }

            ancestors[node_idx.index()] = node_ancestors;
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
        let mut s_values = FxHashMap::with_capacity_and_hasher(dag_nodes.len(), FxBuildHasher);

        s_values.insert(term_idx, 1.0);

        // Process in topological order (children before parents)
        for &node_idx in &self.topo_order {
            if !dag_nodes.contains(&node_idx) || node_idx == term_idx {
                continue;
            }

            let mut max_contribution: f64 = 0.0;
            for child_idx in self.graph.neighbors_directed(node_idx, petgraph::Outgoing) {
                if dag_nodes.contains(&child_idx) {
                    if let Some(&child_s_value) = s_values.get(&child_idx) {
                        max_contribution = max_contribution.max(w * child_s_value);
                    }
                }
            }

            if max_contribution > 0.0 {
                s_values.insert(node_idx, max_contribution);
            }
        }

        s_values
    }

    // /// Calculate Wang similarity between two terms with caching
    // pub fn wang_sim(&self, term1: &str, term2: &str, w: f64) -> Option<f64> {
    //     let term_idx_1 = *self.term_to_idx.get(term1)?;
    //     let term_idx_2 = *self.term_to_idx.get(term2)?;

    //     if term_idx_1 == term_idx_2 {
    //         return Some(1.0);
    //     }

    //     self.wang_sim_by_idx(term_idx_1, term_idx_2, w)
    // }

    // /// Calculate Wang similarity between two term indices (internal method)
    // fn wang_sim_by_idx(&self, term_idx_1: NodeIndex, term_idx_2: NodeIndex, w: f64) -> Option<f64> {
    //     let w_key = w.to_bits(); // Convert f64 to u64 for hashing

    //     let s_val_1 = self.get_or_compute_s_values(term_idx_1, w, w_key);
    //     let s_val_2 = self.get_or_compute_s_values(term_idx_2, w, w_key);

    //     let dag1_nodes = &self.ancestors[term_idx_1.index()];
    //     let dag2_nodes = &self.ancestors[term_idx_2.index()];

    //     // Find intersection more efficiently
    //     let (smaller, larger) = if dag1_nodes.len() < dag2_nodes.len() {
    //         (dag1_nodes, dag2_nodes)
    //     } else {
    //         (dag2_nodes, dag1_nodes)
    //     };

    //     let common_nodes: Vec<NodeIndex> = smaller
    //         .iter()
    //         .filter(|node| larger.contains(node))
    //         .cloned()
    //         .collect();

    //     if common_nodes.is_empty() {
    //         return Some(0.0);
    //     }

    //     let sv1: f64 = s_val_1.values().sum();
    //     let sv2: f64 = s_val_2.values().sum();

    //     let numerator: f64 = common_nodes
    //         .iter()
    //         .map(|&node_idx| {
    //             s_val_1.get(&node_idx).unwrap_or(&0.0) + s_val_2.get(&node_idx).unwrap_or(&0.0)
    //         })
    //         .sum();

    //     let denominator = sv1 + sv2;

    //     if denominator > 0.0 {
    //         Some(numerator / denominator)
    //     } else {
    //         Some(0.0)
    //     }
    // }
}

////////////
// Others //
////////////

/// Filter the similarities based on some threshold
pub fn filter_sims_critval<'a>(
    sim_vals: &[f64],
    names: &'a [String],
    threshold: f64,
) -> Vec<OntoSimRes<'a>> {
    let n = names.len();
    let mut results = Vec::new();
    let mut idx = 0;

    for i in 0..n {
        for j in i..n {
            if i != j {
                if idx < sim_vals.len() {
                    let sim = sim_vals[idx];
                    if sim >= threshold {
                        results.push(OntoSimRes {
                            t1: &names[i],
                            t2: &names[j],
                            sim,
                        })
                    }
                }
                idx += 1;
            }
        }
    }

    results
}
