use petgraph::prelude::*;
use petgraph::visit::{IntoEdges, NodeCount, NodeIndexable};
use petgraph::Graph;
use rayon::prelude::*;
use rustc_hash::FxHashMap;

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
    + std::ops::Add<Output = Self>
    + std::ops::Sub<Output = Self>
    + std::ops::Mul<Output = Self>
    + std::ops::Div<Output = Self>
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
    Max,
    Min,
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
    in_edges_offsets: Vec<usize>,
    out_degrees: Vec<f64>,
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
    pub fn from_petgraph<G>(graph: G) -> Self
    where
        G: NodeCount + IntoEdges + NodeIndexable + Sync,
    {
        let node_count = graph.node_count();

        // Build adjacency structure more efficiently
        let mut out_edges: Vec<Vec<usize>> = vec![Vec::new(); node_count];
        let mut in_edges: Vec<Vec<usize>> = vec![Vec::new(); node_count];

        // Single pass through edges and clippy being dumb
        #[allow(clippy::needless_range_loop)]
        for i in 0..node_count {
            let node_id = graph.from_index(i);
            for edge in graph.edges(node_id) {
                let target_idx = graph.to_index(edge.target());
                out_edges[i].push(target_idx);
                in_edges[target_idx].push(i);
            }
        }

        // Flatten in_edges for better cache locality
        let mut in_edges_flat = Vec::new();
        let mut in_edges_offsets = Vec::with_capacity(node_count + 1);
        in_edges_offsets.push(0);

        for node_in_edges in &in_edges {
            in_edges_flat.extend_from_slice(node_in_edges);
            in_edges_offsets.push(in_edges_flat.len());
        }

        let out_degrees: Vec<f64> = out_edges.iter().map(|edges| edges.len() as f64).collect();

        Self {
            node_count,
            in_edges_flat,
            in_edges_offsets,
            out_degrees,
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
        undirected: bool,
    ) -> Self {
        assert_same_len!(from, to);

        let node_count = nodes.len();

        // Create mapping from node names to indices
        let mut name_to_idx = FxHashMap::default();
        for (idx, name) in nodes.iter().enumerate() {
            name_to_idx.insert(name, idx);
        }

        // Build adjacency lists
        let mut out_edges: Vec<Vec<usize>> = vec![Vec::new(); node_count];
        let mut in_edges: Vec<Vec<usize>> = vec![Vec::new(); node_count];

        for (from_name, to_name) in from.iter().zip(to.iter()) {
            let from_idx = *name_to_idx.get(from_name).unwrap();
            let to_idx = *name_to_idx.get(to_name).unwrap();

            out_edges[from_idx].push(to_idx);
            in_edges[to_idx].push(from_idx);

            if undirected {
                out_edges[to_idx].push(from_idx);
                in_edges[from_idx].push(to_idx);
            }
        }

        // Flatten in_edges for better cache locality
        let mut in_edges_flat = Vec::new();
        let mut in_edges_offsets = Vec::with_capacity(node_count + 1);
        in_edges_offsets.push(0);

        for node_in_edges in &in_edges {
            in_edges_flat.extend_from_slice(node_in_edges);
            in_edges_offsets.push(in_edges_flat.len());
        }

        let out_degrees: Vec<f64> = out_edges.iter().map(|edges| edges.len() as f64).collect();

        Self {
            node_count,
            in_edges_flat,
            in_edges_offsets,
            out_degrees,
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
    fn in_edges(&self, node: usize) -> &[usize] {
        let start = self.in_edges_offsets[node];
        let end = self.in_edges_offsets[node + 1];
        &self.in_edges_flat[start..end]
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
///             fast calculations
/// * `damping_factor` - The dampening factor parameter, i.e., the probability
///                      of resetting.
/// * `personalization_vector` - The vector of probabilities for the reset,
///                              making this the personalised page rank.
/// * `nb_iter` - Maximum number of iterations for the personalised page rank.
/// * `tolerance` - Tolerance of the algorithm.
/// * `working_memory` - The `PageRankWorkingMemory` structure to store the old
///                      and new ranks
///
/// ### Returns
///
/// The (normalised) personalised page rank scores.
pub fn personalised_page_rank_optimised(
    graph: &PageRankGraph,
    damping_factor: f64,
    personalization_vector: &[f64],
    nb_iter: usize,
    tolerance: f64,
    working_memory: &mut PageRankWorkingMemory,
) -> Vec<f64> {
    let node_count = graph.node_count;

    // Reuse pre-allocated vectors
    working_memory.ensure_capacity(node_count);
    let ranks = &mut working_memory.ranks;
    let new_ranks = &mut working_memory.new_ranks;

    // Initialize ranks
    ranks[..node_count].copy_from_slice(personalization_vector);

    let teleport_factor = 1.0 - damping_factor;

    for _ in 0..nb_iter {
        // Compute new ranks
        new_ranks[..node_count]
            .par_iter_mut()
            .enumerate()
            .for_each(|(v, new_rank)| {
                let teleport_prob = teleport_factor * personalization_vector[v];

                let link_prob: f64 = graph
                    .in_edges(v)
                    .iter()
                    .map(|&w| {
                        if graph.out_degrees[w] > 0.0 {
                            damping_factor * ranks[w] / graph.out_degrees[w]
                        } else {
                            damping_factor * ranks[w] * personalization_vector[v]
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
        // This is a crazy command and still frightens me
        std::mem::swap(ranks, new_ranks);

        if squared_norm_2 <= tolerance {
            break;
        }
    }

    // Normalize (make sure tha sum == 1)
    let sum: f64 = ranks[..node_count].iter().sum();
    if sum > 0.0 {
        ranks[..node_count].iter_mut().for_each(|x| *x /= sum);
    }

    ranks[..node_count].to_vec()
}

/// Parallel Personalized Page Rank algorithm.
///
/// ### Params
///
/// * `graph` - The PetGraph on which to run the personalised page-rank.
/// * `damping_factor` - The dampening factor parameter, i.e., the probability
///                      of resetting.
/// * `personalization_vector` - The vector of probabilities for the reset,
///                              making this the personalised page rank.
/// * `nb_iter` - Maximum number of iterations for the personalised page rank.
/// * `tolerance` - Optional tolerance for the algorithm. If not provided, it will
///                 default to `1e-6`.
///
/// ### Returns
///
/// The (normalised) personalised page rank scores.
pub fn personalized_page_rank<G, D>(
    graph: G,
    damping_factor: D,
    personalization_vector: &[D],
    nb_iter: usize,
    tol: Option<D>,
) -> Vec<D>
where
    G: NodeCount + IntoEdges + NodeIndexable + Sync,
    D: NumericType + std::iter::Sum + std::ops::DivAssign,
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
        personalization_vector.len(),
        node_count,
        "Personalization vector length must match node count."
    );

    let tolerance = tol.unwrap_or(D::default_tolerance());

    let mut out_edges: Vec<Vec<usize>> = vec![Vec::new(); node_count];
    let mut in_edges: Vec<Vec<usize>> = vec![Vec::new(); node_count];

    for (i, out_edge_vec) in out_edges.iter_mut().enumerate().take(node_count) {
        let node_id = graph.from_index(i);
        for edge in graph.edges(node_id) {
            let target_idx = graph.to_index(edge.target());
            out_edge_vec.push(target_idx);
            in_edges[target_idx].push(i);
        }
    }

    let out_degrees: Vec<D> = out_edges
        .iter()
        .map(|edges| D::from_f64(edges.len() as f64).unwrap_or(D::zero()))
        .collect();

    let mut ranks: Vec<D> = personalization_vector.to_vec();
    let teleport_factor = D::one() - damping_factor;

    for _ in 0..nb_iter {
        let new_ranks: Vec<D> = (0..node_count)
            .into_par_iter()
            .map(|v| {
                let teleport_prob = teleport_factor * personalization_vector[v];

                let link_prob = in_edges[v]
                    .iter()
                    .map(|&w| {
                        if out_degrees[w] > D::zero() {
                            damping_factor * ranks[w] / out_degrees[w]
                        } else {
                            damping_factor * ranks[w] * personalization_vector[v]
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

/////////////
// Helpers //
/////////////

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
pub fn graph_from_strings(
    nodes: &[String],
    from: &[String],
    to: &[String],
    undirected: bool,
) -> Graph<String, ()> {
    assert_same_len!(from, to);

    let mut graph: Graph<String, ()> = Graph::new();

    let mut term_to_idx = FxHashMap::default();

    for term in nodes {
        let idx = graph.add_node(term.clone());
        term_to_idx.insert(term, idx);
    }

    for (from, to) in from.iter().zip(to.iter()) {
        let from = *term_to_idx.get(from).unwrap();
        let to = *term_to_idx.get(to).unwrap();
        graph.add_edge(from, to, ());
        if undirected {
            graph.add_edge(to, from, ());
        }
    }

    graph
}
