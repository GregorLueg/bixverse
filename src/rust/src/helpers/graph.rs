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

////////////////////////////
// Personalised page rank //
////////////////////////////

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
    D: NumericType + std::iter::Sum,
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

        if squared_norm_2 <= tolerance {
            return new_ranks;
        }

        ranks = new_ranks;
    }

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
