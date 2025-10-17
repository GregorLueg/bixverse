use petgraph::prelude::*;
use petgraph::visit::NodeIndexable;
use petgraph::Graph;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

use std::iter::Sum;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, Sub};

use crate::assert_same_len;

////////////
// Traits //
////////////

/// Simple numeric trait for PageRank calculations
///
/// Needed because of all the trait craziness going on
pub trait NumericType:
    Copy
    + Send
    + Sync
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + PartialOrd
    + std::fmt::Debug
{
    fn zero() -> Self;
    fn one() -> Self;
    fn from_f64(n: f64) -> Option<Self>;
    fn default_tolerance() -> Self;
}

/// Implementation `NumericType` for `f64`
impl NumericType for f64 {
    fn zero() -> Self {
        0.0
    }
    fn one() -> Self {
        1.0
    }
    fn from_f64(n: f64) -> Option<Self> {
        Some(n)
    }
    fn default_tolerance() -> Self {
        1e-6
    }
}

// Implementation `NumericType` for `f64`
impl NumericType for f32 {
    fn zero() -> Self {
        0.0
    }
    fn one() -> Self {
        1.0
    }
    fn from_f64(n: f64) -> Option<Self> {
        Some(n as f32)
    }
    fn default_tolerance() -> Self {
        1e-6
    }
}

///////////
// Enums //
///////////

/// Enum for the TiedSumType types
#[derive(Clone, Debug)]
pub enum TiedSumType {
    /// Maximum between the two diffusion vectors
    Max,
    /// Minimum between the two diffusion vectors
    Min,
    /// Average between the two diffusion vectors
    Avg,
}

/// Parsing the tied summarisation types
///
/// ### Params
///
/// * `s` - The string that defines the tied summarisation type
///
/// ### Results
///
/// The `TiedSumType` defining the tied summarisation type
pub fn parse_tied_sum(s: &str) -> Option<TiedSumType> {
    match s.to_lowercase().as_str() {
        "max" => Some(TiedSumType::Max),
        "min" => Some(TiedSumType::Min),
        "mean" => Some(TiedSumType::Avg),
        _ => None,
    }
}

////////////////////////////
// Personalised page rank //
////////////////////////////

/// Structure for Page Rank Memory
///
/// Allows for faster, better usage of memory
///
/// ### Fields
///
/// * `ranks` - The old ranks
/// * `new_ranks` - The new ranks
#[derive(Debug)]
pub struct PageRankWorkingMemory {
    ranks: Vec<f64>,
    new_ranks: Vec<f64>,
}

impl PageRankWorkingMemory {
    /// Initialise the structure
    ///
    /// ### Returns
    ///
    /// Initialised `PageRankWorkingMemory` structure
    pub fn new() -> Self {
        Self {
            ranks: Vec::new(),
            new_ranks: Vec::new(),
        }
    }

    /// Ensure that the capacity is correct to avoid panics
    ///
    /// ### Params
    ///
    /// * `node_count` - The new node count.
    fn ensure_capacity(&mut self, node_count: usize) {
        if self.ranks.len() < node_count {
            self.ranks.resize(node_count, 0.0);
            self.new_ranks.resize(node_count, 0.0);
        }
    }
}

/// Precomputed graph structure for efficient PageRank computation
///
/// ### Fields
///
/// * `node_count` - The total number of nodes in the graph
/// * `in_edges_flat` - Flattened incoming edges: `[node0_in_edges..., node1_in_edges..., ...]`
/// * `in_edges_offsets` - Offsets into the `in_edges_flat` for each node
/// * `out_degrees` - The out degree for each node.
#[derive(Clone)]
pub struct PageRankGraph {
    node_count: usize,
    in_edges_flat: Vec<usize>,
    in_edge_weights_flat: Vec<f64>,
    in_edges_offsets: Vec<usize>,
    out_weight_sums: Vec<f64>,
}

#[allow(dead_code)]
impl PageRankGraph {
    /// Generate the structure from a given petgraph.
    ///
    /// ### Params
    ///
    /// * `graph` The PetGraph from which to generate the structure.
    ///
    /// ### Returns
    ///
    /// Initialised `PageRankGraph` structure
    pub fn from_petgraph(graph: Graph<&str, f64>) -> Self {
        let node_count = graph.node_count();

        // Build adjacency structure more efficiently
        let mut out_edges: Vec<Vec<(usize, f64)>> = vec![Vec::new(); node_count];
        let mut in_edges: Vec<Vec<(usize, f64)>> = vec![Vec::new(); node_count];

        // Single pass through edges and clippy being dumb

        for i in 0..node_count {
            let node_id = graph.from_index(i);
            for edge in graph.edges(node_id) {
                let target_idx = graph.to_index(edge.target());
                let weight = *edge.weight();
                out_edges[i].push((target_idx, weight));
                in_edges[target_idx].push((i, weight));
            }
        }

        // Flatten in_edges for better cache locality
        let mut in_edges_flat = Vec::new();
        let mut in_edge_weights_flat = Vec::new();
        let mut in_edges_offsets = Vec::with_capacity(node_count + 1);
        in_edges_offsets.push(0);

        for node_in_edges in &in_edges {
            for &(node_idx, weight) in node_in_edges {
                in_edges_flat.push(node_idx);
                in_edge_weights_flat.push(weight);
            }
            in_edges_offsets.push(in_edges_flat.len());
        }

        // Calculate sum of outgoing weights for each node
        let out_weight_sums: Vec<f64> = out_edges
            .iter()
            .map(|edges| edges.iter().map(|(_, weight)| weight).sum())
            .collect();

        Self {
            node_count,
            in_edges_flat,
            in_edge_weights_flat,
            in_edges_offsets,
            out_weight_sums,
        }
    }

    /// Generate the structure directly from node names and edge lists
    ///
    /// ### Params
    ///
    /// * `nodes` - Slice of the node names
    /// * `from` - Slice of the names of the from nodes
    /// * `to` - Slice of the names of the to nodes
    /// * `undirected` - Whether to create bidirectional edges
    ///
    /// ### Returns
    ///
    /// Initialised `PageRankGraph` structure
    pub fn from_strings(
        nodes: &[String],
        from: &[String],
        to: &[String],
        weights: Option<&[f64]>,
        undirected: bool,
    ) -> Self {
        assert_same_len!(from, to);

        if let Some(weights) = weights {
            assert_same_len!(from, weights);
        }

        let node_count = nodes.len();

        // Create mapping from node names to indices
        let mut name_to_idx = FxHashMap::default();
        for (idx, name) in nodes.iter().enumerate() {
            name_to_idx.insert(name, idx);
        }

        // Build adjacency lists with weights
        let mut out_edges: Vec<Vec<(usize, f64)>> = vec![Vec::new(); node_count];
        let mut in_edges: Vec<Vec<(usize, f64)>> = vec![Vec::new(); node_count];

        for (i, (from_name, to_name)) in from.iter().zip(to.iter()).enumerate() {
            let from_idx = *name_to_idx.get(from_name).unwrap();
            let to_idx = *name_to_idx.get(to_name).unwrap();

            let weight = weights.map(|w| w[i]).unwrap_or(1.0);

            out_edges[from_idx].push((to_idx, weight));
            in_edges[to_idx].push((from_idx, weight));

            if undirected {
                out_edges[to_idx].push((from_idx, weight));
                in_edges[from_idx].push((to_idx, weight));
            }
        }

        // Flatten in_edges for better cache locality
        let mut in_edges_flat = Vec::new();
        let mut in_edge_weights_flat = Vec::new();
        let mut in_edges_offsets = Vec::with_capacity(node_count + 1);
        in_edges_offsets.push(0);

        for node_in_edges in &in_edges {
            for &(node_idx, weight) in node_in_edges {
                in_edges_flat.push(node_idx);
                in_edge_weights_flat.push(weight);
            }
            in_edges_offsets.push(in_edges_flat.len());
        }

        // Calculate sum of outgoing weights for each node
        let out_weight_sums: Vec<f64> = out_edges
            .iter()
            .map(|edges| edges.iter().map(|(_, weight)| weight).sum())
            .collect();

        Self {
            node_count,
            in_edges_flat,
            in_edge_weights_flat,
            in_edges_offsets,
            out_weight_sums,
        }
    }

    /// Get incoming edges for a node
    ///
    /// Inline function to hopefully optimise further the compilation of the
    /// program
    ///
    /// ### Params
    ///
    /// * `node` - Get the in_edges for a given node index.
    ///
    /// ### Return
    ///
    /// Returns a slice of in_edges
    #[inline]
    fn in_edges(&self, node: usize) -> (&[usize], &[f64]) {
        let start = self.in_edges_offsets[node];
        let end = self.in_edges_offsets[node + 1];
        (
            &self.in_edges_flat[start..end],
            &self.in_edge_weights_flat[start..end],
        )
    }
}

/// Optimised PageRank with pre-allocated working memory
///
/// This is a highly optimised version of the personalised page rank to be used
/// for rapid permutations. Otherwise, you can just use the other version
///
/// ### Params
///
/// * `graph` - The `PageRankGraph` structure with pre-computed values for
///   fast calculations
/// * `damping_factor` - The dampening factor parameter, i.e., the probability
///   of resetting.
/// * `personalization_vector` - The vector of probabilities for the reset,
///   making this the personalised page rank.
/// * `nb_iter` - Maximum number of iterations for the personalised page rank.
/// * `tolerance` - Tolerance of the algorithm.
/// * `working_memory` - The `PageRankWorkingMemory` structure to store the old
///   and new ranks
///
/// ### Returns
///
/// The (normalised) personalised page rank scores.
pub fn personalised_page_rank_optimised(
    graph: &PageRankGraph,
    damping_factor: f64,
    personalisation_vector: &[f64],
    nb_iter: usize,
    tolerance: f64,
    working_memory: &mut PageRankWorkingMemory,
) -> Vec<f64> {
    let node_count = graph.node_count;

    // reuse pre-allocated vectors
    working_memory.ensure_capacity(node_count);
    let ranks = &mut working_memory.ranks;
    let new_ranks = &mut working_memory.new_ranks;

    // initialise ranks
    ranks[..node_count].copy_from_slice(personalisation_vector);

    let teleport_factor = 1.0 - damping_factor;

    for _ in 0..nb_iter {
        // compute new ranks
        new_ranks[..node_count]
            .par_iter_mut()
            .enumerate()
            .for_each(|(v, new_rank)| {
                let teleport_prob = teleport_factor * personalisation_vector[v];

                let (in_nodes, in_weights) = graph.in_edges(v);
                let link_prob: f64 = in_nodes
                    .iter()
                    .zip(in_weights.iter())
                    .map(|(&w, &edge_weight)| {
                        if graph.out_weight_sums[w] > 0.0 {
                            damping_factor * ranks[w] * edge_weight / graph.out_weight_sums[w]
                        } else {
                            damping_factor * ranks[w] * personalisation_vector[v]
                        }
                    })
                    .sum();

                *new_rank = teleport_prob + link_prob;
            });

        // Check convergence
        let squared_norm_2: f64 = new_ranks[..node_count]
            .par_iter()
            .zip(&ranks[..node_count])
            .map(|(new, old)| {
                let diff = *new - *old;
                diff * diff
            })
            .sum();

        // Swap vectors (no allocation)
        std::mem::swap(ranks, new_ranks);

        if squared_norm_2 <= tolerance {
            break;
        }
    }

    // Normalize (make sure that sum == 1)
    let sum: f64 = ranks[..node_count].iter().sum();
    if sum > 0.0 {
        ranks[..node_count].iter_mut().for_each(|x| *x /= sum);
    }

    ranks[..node_count].to_vec()
}

/// Parallel personalised PageRank algorithm.
///
/// ### Params
///
/// * `graph` - The PetGraph on which to run the personalised page-rank.
/// * `damping_factor` - The dampening factor parameter, i.e., the probability
///   of resetting.
/// * `personalization_vector` - The vector of probabilities for the reset,
///   making this the personalised page rank.
/// * `nb_iter` - Maximum number of iterations for the personalised page rank.
/// * `tolerance` - Optional tolerance for the algorithm. If not provided, it will
///   default to `1e-6`.
///
/// ### Returns
///
/// The (normalised) personalised page rank scores.
pub fn personalised_page_rank<D>(
    graph: Graph<&str, D>,
    damping_factor: D,
    personalisation_vector: &[D],
    nb_iter: usize,
    tol: Option<D>,
) -> Vec<D>
where
    D: NumericType + Sum + DivAssign,
{
    let node_count = graph.node_count();
    if node_count == 0 {
        return vec![];
    }

    // Validate inputs (same as before)
    assert!(
        D::zero() <= damping_factor && damping_factor <= D::one(),
        "Damping factor should be between 0 and 1."
    );
    assert_eq!(
        personalisation_vector.len(),
        node_count,
        "Personalisation vector length must match node count."
    );

    let tolerance = tol.unwrap_or(D::default_tolerance());

    let mut out_edges: Vec<Vec<(usize, D)>> = vec![Vec::new(); node_count];
    let mut in_edges: Vec<Vec<(usize, D)>> = vec![Vec::new(); node_count];

    // build adjacency lists with weights
    for (i, out_edge_vec) in out_edges.iter_mut().enumerate().take(node_count) {
        let node_id = graph.from_index(i);
        for edge in graph.edges(node_id) {
            let target_idx = graph.to_index(edge.target());
            let weight = *edge.weight();
            out_edge_vec.push((target_idx, weight));
            in_edges[target_idx].push((i, weight));
        }
    }

    let out_weight_sums: Vec<D> = out_edges
        .iter()
        .map(|edges| edges.iter().map(|(_, weight)| *weight).sum())
        .collect();

    let mut ranks: Vec<D> = personalisation_vector.to_vec();
    let teleport_factor = D::one() - damping_factor;

    for _ in 0..nb_iter {
        let new_ranks: Vec<D> = (0..node_count)
            .into_par_iter()
            .map(|v| {
                let teleport_prob = teleport_factor * personalisation_vector[v];

                let link_prob = in_edges[v]
                    .iter()
                    .map(|&(w, edge_weight)| {
                        if out_weight_sums[w] > D::zero() {
                            damping_factor * ranks[w] * edge_weight / out_weight_sums[w]
                        } else {
                            damping_factor * ranks[w] * personalisation_vector[v]
                        }
                    })
                    .sum::<D>();

                teleport_prob + link_prob
            })
            .collect();

        let squared_norm_2 = new_ranks
            .par_iter()
            .zip(&ranks)
            .map(|(new, old)| (*new - *old) * (*new - *old))
            .sum::<D>();

        ranks = new_ranks;

        if squared_norm_2 <= tolerance {
            break;
        }
    }

    let sum: D = ranks.iter().copied().sum();
    ranks.iter_mut().for_each(|x| *x /= sum);

    ranks
}

////////////////////////////////////////
// Constrained personalised page rank //
////////////////////////////////////////

/// Constrained parallel personalised PageRank algorithm
///
/// ### Params
///
/// * `graph` - The PetGraph with NodeData and EdgeData
/// * `damping_factor` - The dampening factor parameter, i.e., the probability
///   of resetting.
/// * `personalization_vector` - The vector of probabilities for the reset,
///   making this the personalised page rank.
/// * `nb_iter` - Maximum number of iterations for the personalised page rank.
/// * `tolerance` - Optional tolerance for the algorithm. If not provided, it
///   will default to `1e-6`.
/// * `sink_node_types` - Optional HashSet of node types that act as sinks
///   (force reset of the surfer)
/// * `constrained_edge_types` - Optional HashSet of edge types that force reset
///   after traversal of that edge.
///     
///
/// ### Returns
///
/// The normalised personalised PageRank scores
pub fn constrained_personalised_page_rank<D>(
    graph: &Graph<NodeData, EdgeData>,
    damping_factor: D,
    personalisation_vector: &[D],
    nb_iter: usize,
    tol: Option<D>,
    sink_node_types: Option<&FxHashSet<String>>,
    constrained_edge_types: Option<&FxHashSet<String>>,
) -> Vec<D>
where
    D: NumericType + Sum + Send + Sync + DivAssign + AddAssign,
{
    let node_count = graph.node_count();
    if node_count == 0 {
        return vec![];
    }

    // further assertions
    assert!(
        D::zero() <= damping_factor && damping_factor <= D::one(),
        "Damping factor should be between 0 and 1."
    );
    assert_eq!(
        personalisation_vector.len(),
        node_count,
        "Personalization vector length must match node count."
    );

    let tolerance = tol.unwrap_or(D::default_tolerance());
    let binding = FxHashSet::default();
    let sink_types = sink_node_types.unwrap_or(&binding);
    let constrained_types = constrained_edge_types.unwrap_or(&binding);

    // build transition structure with constraints and weights
    let mut out_edges: Vec<Vec<(usize, f64, bool)>> = vec![Vec::new(); node_count];
    let mut in_edges: Vec<Vec<(usize, f64, bool)>> = vec![Vec::new(); node_count];

    for node_idx in graph.node_indices() {
        let node_idx_usize = graph.to_index(node_idx);
        let node_data = &graph[node_idx];

        // check if this node is a sink - if so, it has no valid outgoing edges
        if sink_types.contains(node_data.node_type) {
            continue;
        }

        // add all outgoing edges with their weights and sink edge flags
        for edge in graph.edges(node_idx) {
            let edge_data = edge.weight();
            let target_idx = graph.to_index(edge.target());
            let is_sink_edge = constrained_types.contains(edge_data.edge_type);

            out_edges[node_idx_usize].push((target_idx, *edge_data.weight, is_sink_edge));
            in_edges[target_idx].push((node_idx_usize, *edge_data.weight, is_sink_edge));
        }
    }

    // calculate out-degrees (sum of weights of ALL edges, including sink edges)
    let out_degrees: Vec<D> = out_edges
        .iter()
        .map(|edges| {
            let total_weight: f64 = edges.iter().map(|(_, weight, _)| weight).sum();
            D::from_f64(total_weight).unwrap_or(D::zero())
        })
        .collect();

    // track mass that can flow forward (non-sink mass) vs absorbed mass (sink mass)
    let mut flowable_ranks: Vec<D> = personalisation_vector.to_vec();
    let mut absorbed_ranks: Vec<D> = vec![D::zero(); node_count];
    let teleport_factor = D::one() - damping_factor;

    for _ in 0..nb_iter {
        // Calculate new flowable and absorbed mass
        let (new_flowable, new_absorbed): (Vec<D>, Vec<D>) = (0..node_count)
            .into_par_iter()
            .map(|v| {
                let teleport_mass = teleport_factor * personalisation_vector[v];
                let mut flowable = teleport_mass;
                let mut absorbed = D::zero();

                for &(w, weight, is_sink_edge) in &in_edges[v] {
                    if out_degrees[w] > D::zero() {
                        let edge_weight = D::from_f64(weight).unwrap_or(D::zero());
                        let flow =
                            damping_factor * flowable_ranks[w] * edge_weight / out_degrees[w];

                        if is_sink_edge {
                            absorbed += flow;
                        } else {
                            flowable += flow;
                        }
                    } else {
                        // Source node w has no outgoing edges (sink node)
                        flowable += damping_factor * flowable_ranks[w] * personalisation_vector[v];
                    }
                }

                (flowable, absorbed)
            })
            .unzip();

        // Total absorbed mass gets teleported back
        let total_absorbed_mass: D = new_absorbed.iter().copied().sum();

        let final_flowable: Vec<D> = new_flowable
            .into_iter()
            .enumerate()
            .map(|(v, flowable)| flowable + total_absorbed_mass * personalisation_vector[v])
            .collect();

        // Total ranks = flowable + absorbed
        let total_ranks: Vec<D> = final_flowable
            .iter()
            .zip(&new_absorbed)
            .map(|(f, a)| *f + *a)
            .collect();

        // Check for convergence using total ranks
        let squared_norm_2 = total_ranks
            .par_iter()
            .zip(flowable_ranks.par_iter().zip(&absorbed_ranks))
            .map(|(new_total, (old_flow, old_abs))| {
                let old_total = *old_flow + *old_abs;
                (*new_total - old_total) * (*new_total - old_total)
            })
            .sum::<D>();

        flowable_ranks = final_flowable;
        absorbed_ranks = new_absorbed;

        if squared_norm_2 <= tolerance {
            break;
        }
    }

    // final ranks = flowable + absorbed
    let mut ranks: Vec<D> = flowable_ranks
        .iter()
        .zip(&absorbed_ranks)
        .map(|(f, a)| *f + *a)
        .collect();

    // normalise
    let sum: D = ranks.iter().copied().sum();
    if sum > D::zero() {
        ranks.iter_mut().for_each(|x| *x /= sum);
    }

    ranks
}

/////////////////////////
// PetGraph generation //
/////////////////////////

/// NodeData structure
///
/// ### Fields
///
/// * `name` - name of the node.
/// * `node_type` - type of the node
#[derive(Debug, Clone)]
#[allow(dead_code)] // clippy is wrongly complaining here
pub struct NodeData<'a> {
    pub name: &'a str,
    pub node_type: &'a str,
}

/// EdgeData structure
///
/// ### Fields
///
/// * `edge_type` - type of the edge.
/// * `weight` - weight of the edge.
#[derive(Debug, Clone)]
pub struct EdgeData<'a> {
    pub edge_type: &'a str,
    pub weight: &'a f64,
}

/// Generate a weighted, labelled PetGraph Graph
///
/// ### Params
///
/// * `nodes` - Slice of node names.
/// * `node_types` - Slice of the node types.
/// * `from` - Slice of the names of the from nodes.
/// * `to` - Slice of the names of the no nodes.
/// * `edge_types` - Slice of the edge types.
/// * `edge_weights` - Slice of the edge weights.
///
/// ### Returns
///
/// A `Graph<NodeData<'_>, EdgeData<'_>>` graph for subsequent usage in
/// constraint personalised page-rank iteration.
///
/// ### Panics
///
/// Function will panic if `nodes` and `nodes_types` do not have the same
/// length and/or when `from`, `to`, `edge_types` and `edge_weights` do not
/// have the same length.
pub fn graph_from_strings_with_attributes<'a>(
    nodes: &'a [String],
    node_types: &'a [String],
    from: &'a [String],
    to: &'a [String],
    edge_types: &'a [String],
    edge_weights: &'a [f64],
) -> Graph<NodeData<'a>, EdgeData<'a>> {
    assert_same_len!(nodes, node_types);
    assert_same_len!(from, to, edge_types, edge_weights);

    let mut graph: Graph<NodeData<'_>, EdgeData<'_>> = Graph::new();
    let mut name_to_idx = FxHashMap::default();

    // add the nodes
    for (name, node_type) in nodes.iter().zip(node_types.iter()) {
        let node_data = NodeData { name, node_type };
        let idx = graph.add_node(node_data);
        name_to_idx.insert(name, idx);
    }

    // add the edges
    for (((from_name, to_name), edge_type), weight) in from
        .iter()
        .zip(to.iter())
        .zip(edge_types.iter())
        .zip(edge_weights.iter())
    {
        let from_idx = *name_to_idx
            .get(from_name)
            .unwrap_or_else(|| panic!("From node '{}' not found in nodes list", from_name));

        let to_idx = *name_to_idx
            .get(to_name)
            .unwrap_or_else(|| panic!("From node '{}' not found in nodes list", from_name));

        let edge_data = EdgeData { edge_type, weight };

        graph.add_edge(from_idx, to_idx, edge_data);
    }

    graph
}

/// Generate a PetGraph Graph
///
/// ### Params
///
/// * `nodes` - Slice of the node names.
/// * `from` - Slice of the names of the from nodes.
/// * `to` - Slice of the names of the no nodes.
/// * `undirected` - Shall a directed or undirected graph be generated.
///
/// ### Returns
///
/// The generated PetGraph
pub fn graph_from_strings<'a>(
    nodes: &'a [String],
    from: &[String],
    to: &[String],
    weights: Option<&[f64]>,
    undirected: bool,
) -> Graph<&'a str, f64> {
    assert_same_len!(from, to);

    if let Some(weights) = weights {
        assert_same_len!(from, weights);
    }

    let mut graph: Graph<&'a str, f64> = Graph::new();

    let mut term_to_idx = FxHashMap::default();

    for term in nodes {
        let idx = graph.add_node(term);
        term_to_idx.insert(term, idx);
    }

    for (i, (from, to)) in from.iter().zip(to.iter()).enumerate() {
        let from = *term_to_idx.get(from).unwrap();
        let to = *term_to_idx.get(to).unwrap();
        let weight = weights.map(|w| w[i]).unwrap_or(1.0);
        graph.add_edge(from, to, weight);
        if undirected {
            graph.add_edge(to, from, weight);
        }
    }

    graph
}

///////////
// Tests //
///////////

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constrained_pagerank() {
        // create a simple test graph: A -> B -> C (where C is a sink node)
        let nodes = vec!["A".to_string(), "B".to_string(), "C".to_string()];
        let node_types = vec![
            "regular".to_string(),
            "regular".to_string(),
            "sink".to_string(),
        ];
        let from = vec!["A".to_string(), "B".to_string()];
        let to = vec!["B".to_string(), "C".to_string()];
        let edge_types = vec!["normal".to_string(), "normal".to_string()];
        let edge_weights = vec![1.0, 1.0];

        let graph = graph_from_strings_with_attributes(
            &nodes,
            &node_types,
            &from,
            &to,
            &edge_types,
            &edge_weights,
        );

        let personalisation = vec![1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];
        let mut sink_types = FxHashSet::default();
        sink_types.insert("sink".to_string());

        let ranks = constrained_personalised_page_rank(
            &graph,
            0.85,
            &personalisation,
            100,
            None,
            Some(&sink_types),
            None,
        );

        println!("Ranks: {:?}", ranks);

        // C (sink node) should have higher rank than A because mass flows A->B->C
        // and C cannot pass mass forward
        assert!(ranks[2] > ranks[0]);

        // Verify the ranks sum to 1
        let sum: f64 = ranks.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_constrained_edge_types() {
        // Test edge type constraints: A -> B -[sink_edge]-> C
        let nodes = vec!["A".to_string(), "B".to_string(), "C".to_string()];
        let node_types = vec!["regular".to_string(); 3];
        let from = vec!["A".to_string(), "B".to_string()];
        let to = vec!["B".to_string(), "C".to_string()];
        let edge_types = vec!["normal".to_string(), "sink_edge".to_string()];
        let edge_weights = vec![1.0, 1.0];

        let graph = graph_from_strings_with_attributes(
            &nodes,
            &node_types,
            &from,
            &to,
            &edge_types,
            &edge_weights,
        );

        let personalisation = vec![1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];
        let mut blocked_edges = FxHashSet::default();
        blocked_edges.insert("sink_edge".to_string());

        let ranks = constrained_personalised_page_rank(
            &graph,
            0.85,
            &personalisation,
            100,
            None,
            None,
            Some(&blocked_edges),
        );

        println!("Edge-constrained ranks: {:?}", ranks);

        // Mass flowing through sink edge B->C gets teleported back
        // C should have most of the mass. It's the same scenario as above
        assert!(ranks[2] > ranks[0] && ranks[2] > ranks[1]);

        // verify normalisation
        let sum: f64 = ranks.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_mixed_edge_types() {
        // Critical test: Node D can be reached by both sink and non-sink edges
        // A -> B -> D (normal edge)
        // C -> D (sink edge)
        // D -> E (normal edge)

        let nodes = vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "D".to_string(),
            "E".to_string(),
        ];
        let node_types = vec!["regular".to_string(); 5];
        let from = vec![
            "A".to_string(),
            "B".to_string(),
            "C".to_string(),
            "D".to_string(),
        ];
        let to = vec![
            "B".to_string(),
            "D".to_string(),
            "D".to_string(),
            "E".to_string(),
        ];
        let edge_types = vec![
            "normal".to_string(),
            "normal".to_string(),
            "sink_edge".to_string(),
            "normal".to_string(),
        ];
        let edge_weights = vec![1.0, 1.0, 1.0, 1.0];

        let graph = graph_from_strings_with_attributes(
            &nodes,
            &node_types,
            &from,
            &to,
            &edge_types,
            &edge_weights,
        );

        let personalisation = vec![0.2; 5]; // uniform
        let mut blocked_edges = FxHashSet::default();
        blocked_edges.insert("sink_edge".to_string());

        let ranks = constrained_personalised_page_rank(
            &graph,
            0.85,
            &personalisation,
            100,
            None,
            None,
            Some(&blocked_edges),
        );

        println!("Mixed edge types ranks: {:?}", ranks);

        // Mass from B->D (normal) should flow to E
        // Mass from C->D (sink) should teleport back
        // So E should have some mass (reachable via A->B->D->E)
        assert!(ranks[4] > 0.1, "E should be reachable via normal path");

        // D should receive mass from both paths but only forward the normal part
        assert!(
            ranks[3] > 0.1,
            "D should accumulate mass from mixed sources"
        );

        let sum: f64 = ranks.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_edge_weights_affect_transitions() {
        // Test that different edge weights DO change PageRank results
        let nodes = vec!["A".to_string(), "B".to_string(), "C".to_string()];
        let node_types = vec!["regular".to_string(); 3];
        let from = vec!["A".to_string(), "A".to_string()];
        let to = vec!["B".to_string(), "C".to_string()];
        let edge_types = vec!["normal".to_string(); 2];

        // First graph: B edge has 10x weight of C edge
        let b_favored_weights = vec![10.0, 1.0];
        let graph1 = graph_from_strings_with_attributes(
            &nodes,
            &node_types,
            &from,
            &to,
            &edge_types,
            &b_favored_weights,
        );

        // Second graph: C edge has 10x weight of B edge
        let c_favored_weights = vec![1.0, 10.0];
        let graph2 = graph_from_strings_with_attributes(
            &nodes,
            &node_types,
            &from,
            &to,
            &edge_types,
            &c_favored_weights,
        );

        let personalisation = vec![1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];

        let ranks1: Vec<f64> = constrained_personalised_page_rank(
            &graph1,
            0.85,
            &personalisation,
            100,
            None,
            None,
            None,
        );

        let ranks2: Vec<f64> = constrained_personalised_page_rank(
            &graph2,
            0.85,
            &personalisation,
            100,
            None,
            None,
            None,
        );

        println!("B-favored weights ranks: {:?}", ranks1);
        println!("C-favored weights ranks: {:?}", ranks2);

        // Higher weight edges should lead to higher ranks
        assert!(
            ranks1[1] > ranks1[2],
            "B should have higher rank when B edge has higher weight"
        );
        assert!(
            ranks2[2] > ranks2[1],
            "C should have higher rank when C edge has higher weight"
        );

        // Verify normalisation
        for ranks in [&ranks1, &ranks2] {
            let sum: f64 = ranks.iter().sum();
            assert!((sum - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_different_personalisation_vectors() {
        let nodes = vec!["A".to_string(), "B".to_string(), "C".to_string()];
        let node_types = vec!["regular".to_string(); 3];
        let from = vec!["A".to_string(), "B".to_string()];
        let to = vec!["B".to_string(), "C".to_string()];
        let edge_types = vec!["normal".to_string(); 2];
        let edge_weights = vec![1.0, 1.0];

        let graph = graph_from_strings_with_attributes(
            &nodes,
            &node_types,
            &from,
            &to,
            &edge_types,
            &edge_weights,
        );

        let uniform_pers = vec![1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];
        let uniform_ranks =
            constrained_personalised_page_rank(&graph, 0.85, &uniform_pers, 100, None, None, None);

        let a_biased_pers = vec![0.8, 0.1, 0.1];
        let a_biased_ranks =
            constrained_personalised_page_rank(&graph, 0.85, &a_biased_pers, 100, None, None, None);

        let c_biased_pers = vec![0.1, 0.1, 0.8];
        let c_biased_ranks =
            constrained_personalised_page_rank(&graph, 0.85, &c_biased_pers, 100, None, None, None);

        println!("Uniform personalization ranks: {:?}", uniform_ranks);
        println!("A-biased personalization ranks: {:?}", a_biased_ranks);
        println!("C-biased personalization ranks: {:?}", c_biased_ranks);

        // biased personalisation should increase ranks of favored nodes
        assert!(
            a_biased_ranks[0] > uniform_ranks[0],
            "A should benefit from A-biased personalization"
        );
        assert!(
            c_biased_ranks[2] > uniform_ranks[2],
            "C should benefit from C-biased personalization"
        );

        // verify normalisation
        for ranks in [&uniform_ranks, &a_biased_ranks, &c_biased_ranks] {
            let sum: f64 = ranks.iter().sum();
            assert!((sum - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_biological_network_example() {
        let nodes = vec![
            "receptor_a".to_string(),
            "receptor_b".to_string(),
            "kinase_c".to_string(),
            "kinase_d".to_string(),
            "tf_e".to_string(),
            "tf_f".to_string(),
            "target_gene_g".to_string(),
            "target_gene_h".to_string(),
            "target_gene_i".to_string(),
        ];
        let node_types = vec![
            "receptor".to_string(),
            "receptor".to_string(),
            "kinase".to_string(),
            "kinase".to_string(),
            "tf".to_string(),
            "tf".to_string(),
            "target_gene".to_string(),
            "target_gene".to_string(),
            "target_gene".to_string(),
        ];
        let from = vec![
            "receptor_a".to_string(),
            "receptor_a".to_string(),
            "receptor_b".to_string(),
            "kinase_c".to_string(),
            "kinase_d".to_string(),
            "tf_e".to_string(),
            "tf_e".to_string(),
            "tf_f".to_string(),
            "tf_f".to_string(),
            "tf_f".to_string(),
        ];
        let to = vec![
            "kinase_c".to_string(),
            "kinase_d".to_string(),
            "kinase_d".to_string(),
            "tf_e".to_string(),
            "tf_f".to_string(),
            "receptor_a".to_string(),
            "target_gene_g".to_string(),
            "receptor_b".to_string(),
            "target_gene_h".to_string(),
            "target_gene_i".to_string(),
        ];
        let edge_types = vec![
            "activation".to_string(),
            "activation".to_string(),
            "activation".to_string(),
            "phosphorylation".to_string(),
            "phosphorylation".to_string(),
            "tf_activation".to_string(),
            "tf_activation".to_string(),
            "tf_activation".to_string(),
            "tf_activation".to_string(),
            "tf_activation".to_string(),
        ];
        let edge_weights = vec![1.0; 10];

        let graph = graph_from_strings_with_attributes(
            &nodes,
            &node_types,
            &from,
            &to,
            &edge_types,
            &edge_weights,
        );

        let personalisation = vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let mut sink_edges = FxHashSet::default();
        sink_edges.insert("activation".to_string());

        let ranks = constrained_personalised_page_rank(
            &graph,
            0.85,
            &personalisation,
            100,
            None,
            None,
            Some(&sink_edges),
        );

        println!("Biological network ranks: {:?}", ranks);

        // kinase_c and kinase_d should be reachable and have non-zero scores
        assert!(
            ranks[2] > 0.0,
            "kinase_c should be reachable via activation edge"
        );
        assert!(
            ranks[3] > 0.0,
            "kinase_d should be reachable via activation edge"
        );

        // Downstream nodes (tf_e, tf_f, target genes) should have very low scores
        // since activation edges cause immediate reset

        let iter = 4..9;
        for i in iter {
            assert!(
                ranks[i] == 0.0,
                "Downstream nodes should have low scores due to activation edge constraints"
            );
        }

        let sum: f64 = ranks.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);
    }
}
