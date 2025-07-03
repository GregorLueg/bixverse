use petgraph::prelude::*;
use petgraph::visit::{IntoEdges, NodeCount, NodeIndexable};
use petgraph::Graph;
use rayon::prelude::*;
use rustc_hash::FxHashMap;

////////////
// Traits //
////////////

/// Simple numeric trait for PageRank calculations
/// Need this because of all the trait craziness going on
#[allow(dead_code)]
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
    fn from_usize(n: usize) -> Self;
    fn from_f64(n: f64) -> Option<Self>;
    fn abs(self) -> Self;
    fn default_tolerance() -> Self;
}

// Implementation for f64
impl NumericType for f64 {
    fn zero() -> Self {
        0.0
    }
    fn one() -> Self {
        1.0
    }
    fn from_usize(n: usize) -> Self {
        n as f64
    }
    fn from_f64(n: f64) -> Option<Self> {
        Some(n)
    }
    fn abs(self) -> Self {
        self.abs()
    }
    fn default_tolerance() -> Self {
        1e-6
    }
}

// Implementation for f32
impl NumericType for f32 {
    fn zero() -> Self {
        0.0
    }
    fn one() -> Self {
        1.0
    }
    fn from_usize(n: usize) -> Self {
        n as f32
    }
    fn from_f64(n: f64) -> Option<Self> {
        Some(n as f32)
    }
    fn abs(self) -> Self {
        self.abs()
    }
    fn default_tolerance() -> Self {
        1e-6
    }
}

///////////
// Enums //

/// Enum for the ICA types
#[derive(Clone, Debug)]
pub enum TiedSumType {
    Max,
    Min,
    Avg,
}

/// Parsing the tied summarisation types
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

/// Reusable working memory to avoid allocations
/// This is a good suggestion from Claude! Also, usage of memory swap makes stuff
/// way better...
#[derive(Debug)]
pub struct PageRankWorkingMemory {
    ranks: Vec<f64>,
    new_ranks: Vec<f64>,
}

impl PageRankWorkingMemory {
    /// Generaten of a new class
    pub fn new() -> Self {
        Self {
            ranks: Vec::new(),
            new_ranks: Vec::new(),
        }
    }

    /// Ensure that the capacity is correct to avoid panics
    fn ensure_capacity(&mut self, node_count: usize) {
        if self.ranks.len() < node_count {
            self.ranks.resize(node_count, 0.0);
            self.new_ranks.resize(node_count, 0.0);
        }
    }
}

/// Precomputed graph structure for efficient PageRank computation
#[derive(Clone)]
pub struct PageRankGraph {
    node_count: usize,            // Node counts
    in_edges_flat: Vec<usize>, // Flattened incoming edges: [node0_in_edges..., node1_in_edges..., ...]
    in_edges_offsets: Vec<usize>, // Offsets into in_edges_flat for each node
    out_degrees: Vec<f64>,     // Out-degrees for each node
}

impl PageRankGraph {
    /// Generate the structure from a given petgraph.
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

    /// Get incoming edges for a node
    #[inline]
    fn in_edges(&self, node: usize) -> &[usize] {
        let start = self.in_edges_offsets[node];
        let end = self.in_edges_offsets[node + 1];
        &self.in_edges_flat[start..end]
    }
}

/// Optimized PageRank with pre-allocated working memory
pub fn personalized_page_rank_optimized(
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

/// Generate a graph object from three vectors
pub fn graph_from_strings(
    nodes: &[String],
    from: &[String],
    to: &[String],
    undirected: bool,
) -> Graph<String, ()> {
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
